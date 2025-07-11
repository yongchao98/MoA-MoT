def interpret_pcr_result(ct_value, cutoff_cycles, threshold_rfu):
    """
    Interprets a PCR result based on the Ct value and assay cut-off.
    """
    print(f"Analysis of PCR Result:")
    print(f" - Threshold: {threshold_rfu:,} RFU")
    print(f" - Assay Cut-off: {cutoff_cycles} cycles")
    print(f" - Measured Ct Value: {ct_value} cycles")
    print("-" * 30)

    # The core logic is the comparison between the Ct value and the cut-off.
    print(f"To determine the result, we check if the Ct value ({ct_value}) is less than the assay cut-off ({cutoff_cycles}).")

    if ct_value < cutoff_cycles:
        interpretation = "POSITIVE"
        reason = f"The Ct value of {ct_value} is less than the cut-off of {cutoff_cycles}. This indicates that sufficient amplification occurred for the fluorescence signal to cross the threshold."
    else:
        interpretation = "NEGATIVE"
        reason = f"The Ct value ({ct_value}) is not less than the cut-off of {cutoff_cycles} (or no Ct was detected). This indicates the fluorescence signal did not cross the threshold within the specified number of cycles."

    print(f"\nFinal Interpretation: {interpretation}")
    print(f"Reason: {reason}")


# Given parameters from the problem
pcr_threshold = 125000
pcr_cutoff = 40
pcr_ct_value = 30.8

# Run the interpretation
interpret_pcr_result(pcr_ct_value, pcr_cutoff, pcr_threshold)