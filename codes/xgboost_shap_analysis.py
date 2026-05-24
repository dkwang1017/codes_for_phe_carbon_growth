#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
XGBoost-SHAP analysis for drivers of assimilation-growth coupling.

Purpose
-------
This script fits an XGBoost regression model to evaluate the relative importance
and marginal effects of environmental and phenological predictors on the
coupling between gross primary productivity and tree growth (RGPP-Growth).

The workflow includes:
1. Reading the processed site-level input table.
2. Filtering site records and removing invalid or missing predictor values.
3. Splitting the dataset into training and testing subsets.
4. Performing grid-search hyperparameter evaluation using K-fold cross-validation
   within the training subset.
5. Fitting the final XGBoost model using the selected hyperparameters.
6. Evaluating model performance on the held-out test subset.
7. Calculating SHAP values to interpret predictor importance and response patterns.
8. Saving model-performance statistics, SHAP feature importance, selected
   hyperparameters and diagnostic figures.

Notes for editors and reviewers
-------------------------------
This script uses XGBoost as an interpretable machine-learning approach, together
with SHAP, to support ecological attribution of RGPP-Growth variability. It is
not intended to develop an operational forecasting model.

Input
-----
A processed CSV table containing one row per tree-ring site and the following
columns:

Required response variable:
    max_correlation

Required predictors:
    sos_trend, gpp_mean, ssrd_mean, pre_mean, avg_sos, lai_mean, tmp_mean,
    sm_mean, avg_frost_days, sand_mean, oc_mean, clay_mean, AI

Required site identifier:
    site

Output
------
The script writes the following files to output_dir:
    xgboost_model_performance.csv
    xgboost_grid_search_results.csv
    xgboost_selected_hyperparameters.csv
    shap_feature_importance.csv
    shap_summary_plot.png
    shap_feature_importance_barplot.png
    cleaned_model_input_data.csv

Software requirements
---------------------
Python >= 3.9
pandas, numpy, matplotlib, scikit-learn, xgboost, shap
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shap
import xgboost as xgb
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import KFold, ParameterGrid, train_test_split


# =============================================================================
# 1. User-defined paths and parameters
# =============================================================================

# Replace these paths with the local paths to the processed input table and
# desired output folder.
input_file = Path("path/to/site_level_predictors_and_rgpp_growth.csv")
output_dir = Path("path/to/output/xgboost_shap_analysis")
output_dir.mkdir(parents=True, exist_ok=True)

# Response and predictor variables
response_var = "max_correlation"

predictor_vars = [
    "sos_trend",
    "gpp_mean",
    "ssrd_mean",
    "pre_mean",
    "avg_sos",
    "lai_mean",
    "tmp_mean",
    "sm_mean",
    "avg_frost_days",
    "sand_mean",
    "oc_mean",
    "clay_mean",
    "AI",
]

site_col = "site"

# General filtering parameters
invalid_fill_value = 65535
max_aridity_index = 4
max_soil_organic_carbon = 10

# Train-test split and cross-validation settings
test_size = 0.30
random_state = 42
n_splits = 5

# XGBoost hyperparameter grid
param_grid = {
    "max_depth": [3, 5, 7],
    "learning_rate": [0.01, 0.05, 0.10],
    "n_estimators": [500, 1000, 1500],
    "subsample": [0.6, 0.8, 1.0],
    "colsample_bytree": [0.5, 0.7, 1.0],
}

# Number of CPU threads used by XGBoost
n_jobs = 4


# =============================================================================
# 2. Helper functions
# =============================================================================

def filter_valid_site_ids(df: pd.DataFrame, site_column: str) -> pd.DataFrame:
    """Keep standard ITRDB-like site IDs composed of letters followed by digits."""
    if site_column not in df.columns:
        raise ValueError(f"Required site column '{site_column}' was not found.")
    pattern = re.compile(r"^[A-Za-z]+[0-9]+$")
    return df[df[site_column].astype(str).str.match(pattern, na=False)].copy()


def prepare_model_data(df: pd.DataFrame) -> pd.DataFrame:
    """Clean and subset the input table for XGBoost-SHAP analysis."""
    required_cols = [site_col] + predictor_vars + [response_var]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"The input table is missing required columns: {missing_cols}")

    n_initial = len(df)

    # Keep standard site identifiers.
    df = filter_valid_site_ids(df, site_col)
    n_after_site_filter = len(df)

    # Keep only required columns and remove missing predictor/response values.
    data = df[required_cols].copy()
    data = data.dropna(subset=predictor_vars + [response_var])
    n_after_nan_filter = len(data)

    # Remove invalid fill values.
    invalid_mask = data[predictor_vars].eq(invalid_fill_value).any(axis=1)
    data = data.loc[~invalid_mask].copy()
    n_after_fill_filter = len(data)

    # Apply conservative filters to remove implausible auxiliary values.
    data = data[data["AI"] < max_aridity_index].copy()
    data = data[data["oc_mean"] < max_soil_organic_carbon].copy()
    n_final = len(data)

    print("Data-cleaning summary")
    print("---------------------")
    print(f"Initial rows: {n_initial}")
    print(f"Rows after site-ID filtering: {n_after_site_filter}")
    print(f"Rows after removing missing values: {n_after_nan_filter}")
    print(f"Rows after removing invalid fill values: {n_after_fill_filter}")
    print(f"Final rows used for modeling: {n_final}")
    print()

    if n_final < 10:
        raise ValueError("Too few observations remain after filtering.")

    return data


def evaluate_regression(y_true: np.ndarray, y_pred: np.ndarray) -> dict:
    """Return standard regression-performance metrics."""
    return {
        "rmse": float(np.sqrt(mean_squared_error(y_true, y_pred))),
        "r2": float(r2_score(y_true, y_pred)),
        "mae": float(mean_absolute_error(y_true, y_pred)),
    }


def cross_validate_params(
    X_train: pd.DataFrame,
    y_train: pd.Series,
    params: dict,
    n_splits: int,
    random_state: int,
) -> dict:
    """Evaluate one XGBoost hyperparameter set using K-fold cross-validation."""
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=random_state)

    fold_r2 = []
    fold_rmse = []
    fold_mae = []

    for train_idx, val_idx in kf.split(X_train):
        X_tr = X_train.iloc[train_idx]
        X_val = X_train.iloc[val_idx]
        y_tr = y_train.iloc[train_idx]
        y_val = y_train.iloc[val_idx]

        model = xgb.XGBRegressor(
            objective="reg:squarederror",
            random_state=random_state,
            n_jobs=n_jobs,
            **params,
        )
        model.fit(X_tr, y_tr)

        y_val_pred = model.predict(X_val)
        metrics = evaluate_regression(y_val, y_val_pred)

        fold_r2.append(metrics["r2"])
        fold_rmse.append(metrics["rmse"])
        fold_mae.append(metrics["mae"])

    return {
        "mean_cv_r2": float(np.mean(fold_r2)),
        "sd_cv_r2": float(np.std(fold_r2, ddof=1)),
        "mean_cv_rmse": float(np.mean(fold_rmse)),
        "sd_cv_rmse": float(np.std(fold_rmse, ddof=1)),
        "mean_cv_mae": float(np.mean(fold_mae)),
        "sd_cv_mae": float(np.std(fold_mae, ddof=1)),
    }


def select_hyperparameters(X_train: pd.DataFrame, y_train: pd.Series) -> tuple[dict, pd.DataFrame]:
    """
    Select XGBoost hyperparameters using cross-validated predictive performance.

    The selected parameter set maximizes mean cross-validated R². If multiple
    parameter sets have the same R², the one with lower mean RMSE is selected.
    """
    records = []

    for params in ParameterGrid(param_grid):
        cv_metrics = cross_validate_params(
            X_train=X_train,
            y_train=y_train,
            params=params,
            n_splits=n_splits,
            random_state=random_state,
        )
        record = {**params, **cv_metrics}
        records.append(record)

    results = pd.DataFrame(records)
    results = results.sort_values(
        by=["mean_cv_r2", "mean_cv_rmse"],
        ascending=[False, True],
    ).reset_index(drop=True)

    best_params = {
        "max_depth": int(results.loc[0, "max_depth"]),
        "learning_rate": float(results.loc[0, "learning_rate"]),
        "n_estimators": int(results.loc[0, "n_estimators"]),
        "subsample": float(results.loc[0, "subsample"]),
        "colsample_bytree": float(results.loc[0, "colsample_bytree"]),
    }

    return best_params, results


def calculate_shap_outputs(
    model: xgb.XGBRegressor,
    X_for_shap: pd.DataFrame,
    output_dir: Path,
) -> pd.DataFrame:
    """Calculate SHAP values and save feature-importance outputs."""
    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(X_for_shap)

    mean_abs_shap = np.mean(np.abs(shap_values), axis=0)
    shap_importance = pd.DataFrame(
        {
            "feature": X_for_shap.columns,
            "mean_abs_shap": mean_abs_shap,
        }
    ).sort_values("mean_abs_shap", ascending=False)

    shap_importance["rank"] = np.arange(1, len(shap_importance) + 1)
    shap_importance.to_csv(output_dir / "shap_feature_importance.csv", index=False)

    # SHAP summary plot
    plt.figure()
    shap.summary_plot(
        shap_values,
        X_for_shap,
        feature_names=X_for_shap.columns,
        show=False,
    )
    plt.tight_layout()
    plt.savefig(output_dir / "shap_summary_plot.png", dpi=300, bbox_inches="tight")
    plt.close()

    # Mean absolute SHAP feature-importance bar plot
    plot_data = shap_importance.sort_values("mean_abs_shap", ascending=True)
    plt.figure(figsize=(8, 6))
    plt.barh(plot_data["feature"], plot_data["mean_abs_shap"])
    plt.xlabel("Mean absolute SHAP value")
    plt.ylabel("Predictor")
    plt.title("SHAP-based feature importance")
    plt.tight_layout()
    plt.savefig(output_dir / "shap_feature_importance_barplot.png", dpi=300, bbox_inches="tight")
    plt.close()

    return shap_importance


# =============================================================================
# 3. Main workflow
# =============================================================================

def main() -> None:
    """Run the complete XGBoost-SHAP workflow."""
    raw_data = pd.read_csv(input_file)
    model_data = prepare_model_data(raw_data)

    model_data.to_csv(output_dir / "cleaned_model_input_data.csv", index=False)

    X = model_data[predictor_vars].copy()
    y = model_data[response_var].copy()

    # Split the processed dataset into training and held-out testing subsets.
    # No scaling is applied because tree-based models such as XGBoost do not
    # require standardized predictors.
    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=test_size,
        random_state=random_state,
    )

    print(f"Training samples: {len(X_train)}")
    print(f"Testing samples: {len(X_test)}")
    print()

    # Hyperparameter selection by K-fold cross-validation within the training set.
    best_params, grid_results = select_hyperparameters(X_train, y_train)
    grid_results.to_csv(output_dir / "xgboost_grid_search_results.csv", index=False)

    with open(output_dir / "xgboost_selected_hyperparameters.json", "w", encoding="utf-8") as f:
        json.dump(best_params, f, indent=2)

    pd.DataFrame([best_params]).to_csv(
        output_dir / "xgboost_selected_hyperparameters.csv",
        index=False,
    )

    print("Selected hyperparameters:")
    print(best_params)
    print()

    # Fit final model using the selected hyperparameters.
    final_model = xgb.XGBRegressor(
        objective="reg:squarederror",
        random_state=random_state,
        n_jobs=n_jobs,
        **best_params,
    )
    final_model.fit(X_train, y_train)

    # Evaluate the final model on the held-out testing subset.
    y_test_pred = final_model.predict(X_test)
    test_metrics = evaluate_regression(y_test, y_test_pred)

    performance = {
        "n_total": len(model_data),
        "n_train": len(X_train),
        "n_test": len(X_test),
        "test_size": test_size,
        "random_state": random_state,
        **test_metrics,
    }

    pd.DataFrame([performance]).to_csv(
        output_dir / "xgboost_model_performance.csv",
        index=False,
    )

    print("Held-out test performance:")
    for key, value in test_metrics.items():
        print(f"{key}: {value:.4f}")
    print()

    # SHAP interpretation is calculated for the full cleaned dataset to describe
    # predictor importance and marginal response patterns across all available sites.
    shap_importance = calculate_shap_outputs(final_model, X, output_dir)

    print("Top predictors by mean absolute SHAP value:")
    print(shap_importance.head(10).to_string(index=False))
    print()
    print(f"All outputs saved to: {output_dir}")


if __name__ == "__main__":
    main()
