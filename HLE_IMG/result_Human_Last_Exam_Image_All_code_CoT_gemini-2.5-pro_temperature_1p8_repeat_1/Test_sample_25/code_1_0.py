def interpret_pcr_result(ct_value, cutoff_cycle, threshold_rfu):
    """
    Interprets a PCR result based on the Ct value and assay cut-off.

    Args:
        ct_value (float): The cycle at which fluorescence crossed the threshold.
        cutoff_cycle (int): The maximum number of cycles for the assay.
        threshold_rfu (int): The fluorescence threshold in RFU.
    """
    print(f"Interpreting PCR Result:")
    print(f"- Fluorescence Threshold: {threshold_rfu:,} RFU")
    print(f"- Assay Cut-off: {cutoff_cycle} cycles")
    print(f"- Measured Ct Value: {ct_value} cycles")
    print("-" * 25)

    print("Rule: The result is 'Positive' if the Ct value is less than the assay cut-off cycle.")
    print(f"Comparison: Is {ct_value} < {cutoff_cycle}?")

    if ct_value < cutoff_cycle:
        result = "Positive"
        is_less = "Yes"
    else:
        result = "Negative"
        is_less = "No"

    print(f"Result of comparison: {is_less}")
    print("-" * 25)
    print(f"Final Interpretation: The sample is {result}.")
    print("This is because significant amplification was detected before the assay completed its run.")

# Given values from the problem description
pcr_threshold = 125000
pcr_assay_cutoff = 40
pcr_ct_value = 30.8

# Run the interpretation
interpret_pcr_result(pcr_ct_value, pcr_assay_cutoff, pcr_threshold)
