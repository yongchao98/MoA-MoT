def interpret_pcr_result(ct_value, assay_cutoff, threshold_rfu):
    """
    Interprets a PCR result based on the Ct value and assay cut-off.

    Args:
        ct_value (float): The cycle threshold value.
        assay_cutoff (int): The maximum number of cycles for the assay.
        threshold_rfu (int): The fluorescence threshold in RFU.
    """
    print("--- PCR Result Interpretation ---")
    print(f"Threshold (RFU): {threshold_rfu}")
    print(f"Assay Cut-off (Cycles): {assay_cutoff}")
    print(f"Measured Ct Value: {ct_value}")
    print("-" * 33)

    print("The rule for interpretation is to compare the Ct value with the assay cut-off.")
    print(f"Is the Ct value ({ct_value}) less than the assay cut-off ({assay_cutoff})?")

    if ct_value < assay_cutoff:
        print(f"\nYes, {ct_value} is less than {assay_cutoff}.")
        interpretation = "The sample is POSITIVE. The amplification curve crossed the threshold before the assay cut-off."
    else:
        print(f"\nNo, {ct_value} is not less than {assay_cutoff}.")
        interpretation = "The sample is NEGATIVE. The amplification curve did not cross the threshold within the assay's cycle limit."

    print("\nFinal Interpretation:")
    print(interpretation)


# Given values from the problem
pcr_ct_value = 30.8
pcr_assay_cutoff = 40
pcr_threshold_rfu = 125000

# Run the interpretation
interpret_pcr_result(pcr_ct_value, pcr_assay_cutoff, pcr_threshold_rfu)
