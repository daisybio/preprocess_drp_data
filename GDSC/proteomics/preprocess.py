from typing import Union, List

import numpy as np
import pandas as pd


def extract_and_filter_protein_lfq_data(
    df: pd.DataFrame,
    intensity_column_prefix: str = "LFQ intensity ",
    keep_intensity_column_prefix: bool = False,
) -> pd.DataFrame:
    df = df[df["Reverse"] != "+"]
    df = df[df["Potential contaminant"] != "+"]
    if "Only identified by site" in df.columns:
        df = df[df["Only identified by site"] != "+"]

    df = df.set_index("Gene names")
    df = df.filter(like=intensity_column_prefix)
    df.index = df.index.str.split(";").str[0]
    if not keep_intensity_column_prefix:
        df.columns = df.columns.str.removeprefix(intensity_column_prefix)
    return df


def log_transform(df: pd.DataFrame) -> pd.DataFrame:
    df = np.log10(df)
    df = df.replace(-np.inf, np.nan)
    return df


def filter_by_min_occurrence(
    df: pd.DateOffset, min_occurrence_rate: float = 0.7
) -> pd.DataFrame:
    return df[df.count(axis=1) >= min_occurrence_rate * len(df.columns)]


def impute_normal_down_shift_distribution(
    unimputed_df: pd.DataFrame, column_wise=True, width=0.3, downshift=1.8, seed=100
) -> pd.DataFrame:
    """
    Performs imputation across a matrix columnswise
    https://rdrr.io/github/jdreyf/jdcbioinfo/man/impute_normal.html#google_vignette
    :width: Scale factor for the standard deviation of imputed distribution relative to the sample standard deviation.
    :downshift: Down-shifted the mean of imputed distribution from the sample mean, in units of sample standard deviation.
    :seed: Random seed

    """
    unimputed_matrix = unimputed_df.to_numpy()
    columns_names = unimputed_df.columns
    rownames = unimputed_df.index
    unimputed_matrix[~np.isfinite(unimputed_matrix)] = None
    main_mean = np.nanmean(unimputed_matrix)
    main_std = np.nanstd(unimputed_matrix)
    np.random.seed(seed=seed)

    def impute_normal_per_vector(temp: np.ndarray, width=width, downshift=downshift):
        """Performs imputation for a single vector"""
        if column_wise:
            temp_sd = np.nanstd(temp)
            temp_mean = np.nanmean(temp)
        else:
            # over all matrix
            temp_sd = main_std
            temp_mean = main_mean

        shrinked_sd = width * temp_sd
        downshifted_mean = temp_mean - (downshift * temp_sd)
        n_missing = np.count_nonzero(np.isnan(temp))
        temp[np.isnan(temp)] = np.random.normal(
            loc=downshifted_mean, scale=shrinked_sd, size=n_missing
        )
        return temp

    final_matrix = np.apply_along_axis(impute_normal_per_vector, 0, unimputed_matrix)
    final_df = pd.DataFrame(final_matrix)
    final_df.index = rownames
    final_df.columns = columns_names
    return final_df


def median_centering_ms1(merged_ms1_df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalizes samples by multiplying each batch MS1 with its own correction factor.
    Only uses peptides detected in >70% of the batches for computing median to prevent
    low-abundant peptides from dragging down the median in runs with deeper coverage
    :param all_batches: list of evidence dataframes
    :return:
    """
    # the rows are the proteins and columns are samples
    # take the median for each sample after considering the proteins which are at least in more than 70% samples
    medians = (
        merged_ms1_df[merged_ms1_df.count(axis=1) > 0.7 * len(merged_ms1_df.columns)]
        .median(axis=0)
        .to_dict()
    )

    # the average of the medians across all samples
    mean_median = pd.Series(medians.values()).mean()

    df = merged_ms1_df.copy()
    for sample_name in merged_ms1_df.columns:

        # we need a compensation factor to compare different samples with each other
        correction_factor = mean_median / medians[sample_name]
        # print(f"Correction factor for {sample_name}: {round(correction_factor, 3)}")

        df[sample_name] = df[sample_name] * correction_factor
    return df


def normalize_intensity_by_reference(
    row: pd.Series,
    reference_channel: Union[str, int],
    list_of_batches: List[Union[str, int]],
) -> pd.Series:
    """
    Normalizes intensity values in a dataset based on a reference channel for each batch.
    Intensities are scaled back to their original value range using the median of all reference intensity values.

    Parameters:
        row (pd.Series): The data row containing intensity values.
        reference_channel (str): The name of the reference channel to normalize by.
        list_of_batches (list): List of batch identifiers to process.

    Returns:
        pd.Series: The normalized intensity values for the row.
    """

    # Extract intensity values for the reference channel across all batches
    reference_values = row[
        row.index.str.startswith(f"Reporter intensity corrected {reference_channel} ")
    ]

    # Remove zero values and calculate the median of the non-zero reference intensities
    non_zero_references = reference_values[reference_values != 0]
    reference_median = (
        non_zero_references.median() if not non_zero_references.empty else 0
    )

    # Create a mask to identify all channels starting with "Reporter intensity corrected"
    reporter_mask = row.index.str.startswith("Reporter intensity corrected ")

    # Initialize a dictionary to store the normalized values
    res_dict = row.to_dict()

    # Iterate through each batch
    for batch in list_of_batches:
        # Extract reference intensity for the current batch
        reference_col = f"Reporter intensity corrected {reference_channel} {batch}"
        reference_intensity = row.get(reference_col, 0)

        # Identify all channels for the current batch
        batch_channels_mask = reporter_mask & row.index.str.endswith(f" {batch}")
        batch_values = row[batch_channels_mask]

        # Normalize batch values using reference intensity, if non-zero
        if reference_intensity != 0:
            normalized_values = (batch_values / reference_intensity) * reference_median
        else:
            normalized_values = pd.Series(0, index=batch_values.index)

        # Update the result dictionary with the normalized values
        res_dict.update(normalized_values.to_dict())

    return pd.Series(res_dict)
