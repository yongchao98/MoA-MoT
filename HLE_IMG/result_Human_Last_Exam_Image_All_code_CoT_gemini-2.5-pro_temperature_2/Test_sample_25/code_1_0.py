def interpret_pcr_result(ct_value, cutoff_cycles, threshold_rfu):
    """
    Interprets a PCR result based on the Ct value, cycle cutoff, and threshold.

    Args:
        ct_value (float): The cycle threshold value.
        cutoff_cycles (int): The assay's cycle cut-off.
        threshold_rfu (int): The fluorescence threshold in RFU.
    """

    print(f"Given Parameters:")
    print(f"- Fluorescence Threshold: {threshold_rfu:,} RFU")
    print(f"- Assay Cut-off: {cutoff_cycles} cycles")
    print(f"- Observed Ct value: {ct_value}")
    print("\nInterpretation:")

    # Check if the Ct value is below the assay cut-off.
    if ct_value < cutoff_cycles:
        result = "Positive"
        explanation = (
            f"The result is interpreted as '{result}'.\n"
            f"This is because the sample's amplification curve crossed the threshold of {threshold_rfu:,} RFU at cycle {ct_value}, "
            f"which is before the assay cut-off of {cutoff_cycles} cycles.\n"
            "A positive result indicates that the target genetic material was detected in the sample.\n"
            f"The Ct value of {ct_value} is considered relatively high, which suggests a low initial concentration of the target."
        )
    else:
        result = "Negative"
        explanation = (
            f"The result is interpreted as '{result}'.\n"
            f"This is because the sample's amplification curve did not cross the threshold of {threshold_rfu:,} RFU "
            f"before the assay cut-off of {cutoff_cycles} cycles.\n"
            "A negative result indicates that the target genetic material was not detected in the sample."
        )
    
    print(explanation)
    # The final answer is the interpretation text. Since the task asks for a concise answer format,
    # and the primary interpretation is 'Positive', we will use that as the basis.
    # However, the full text explanation is the most helpful output.
    # For a multiple-choice-style final answer, "Positive" would be the key takeaway.


# Parameters from the problem description
threshold = 125000
cutoff = 40
ct = 30.8

interpret_pcr_result(ct, cutoff, threshold)