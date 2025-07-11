def interpret_pcr_result(ct_value, assay_cutoff, threshold):
    """
    Interprets the result of a PCR assay based on Ct value and cut-off.

    Args:
        ct_value (float): The cycle at which fluorescence crossed the threshold.
        assay_cutoff (int): The maximum number of cycles for the assay.
        threshold (int): The fluorescence threshold in RFU.
    """
    print("--- PCR Assay Interpretation ---")
    print(f"Fluorescence Threshold (RFU): {threshold}")
    print(f"Assay Cut-off (Cycles): {assay_cutoff}")
    print(f"Measured Ct Value (Cycles): {ct_value}")
    print("-" * 30)

    # A result is positive if amplification is detected before the cut-off.
    if 0 < ct_value < assay_cutoff:
        result = "Positive"
        explanation = (
            f"The sample is POSITIVE because the Ct value ({ct_value}) is less than the assay cut-off ({assay_cutoff}).\n"
            f"This means the sample's fluorescence crossed the threshold of {threshold} RFU before the end of the assay, "
            "indicating the presence of the target genetic material."
        )
    else:
        result = "Negative"
        explanation = (
            f"The sample is NEGATIVE because the Ct value ({ct_value}) is not within the valid range below the assay cut-off ({assay_cutoff}).\n"
            "This means amplification was not detected in the sample."
        )

    print(explanation)
    return result

# Given parameters from the problem
pcr_threshold = 125000 # As stated in the text, although graph shows 125
pcr_assay_cutoff = 40
pcr_ct_value = 30.8

# Run the interpretation
final_result = interpret_pcr_result(pcr_ct_value, pcr_assay_cutoff, pcr_threshold)
