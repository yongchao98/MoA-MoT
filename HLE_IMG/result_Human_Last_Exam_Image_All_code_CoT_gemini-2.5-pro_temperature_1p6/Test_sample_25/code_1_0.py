def interpret_pcr_result(ct_value, threshold_rfu, cutoff_cycles):
    """
    Interprets the result of a PCR assay based on Ct value, threshold, and cutoff.

    Args:
        ct_value (float): The cycle threshold value.
        threshold_rfu (int): The fluorescence threshold in RFU.
        cutoff_cycles (int): The assay cutoff in cycles.
    """
    print("--- PCR Assay Interpretation ---")
    print(f"Fluorescence Threshold: {threshold_rfu:,} RFU")
    print(f"Assay Cut-off: {cutoff_cycles} cycles")
    print(f"Measured Ct Value: {ct_value}")
    print("---------------------------------")

    # A Ct value is generated only when the curve crosses the threshold.
    # We check if the Ct value is below the assay cutoff.
    if ct_value < cutoff_cycles:
        result = "POSITIVE"
        explanation = (
            f"The amplification curve crossed the threshold of {threshold_rfu:,} RFU "
            f"at cycle {ct_value}. Since this Ct value is less than the assay "
            f"cut-off of {cutoff_cycles} cycles, the target was successfully detected."
        )
    else:
        result = "NEGATIVE"
        explanation = (
            f"The amplification curve did not cross the threshold of {threshold_rfu:,} RFU "
            f"before the assay cut-off of {cutoff_cycles} cycles."
        )

    print(f"Final Interpretation: {result}")
    print(f"Explanation: {explanation}")


# Given parameters from the problem
threshold = 125000
assay_cutoff = 40
ct = 30.8

# Run the interpretation
interpret_pcr_result(ct_value=ct, threshold_rfu=threshold, cutoff_cycles=assay_cutoff)