def interpret_pcr_result(ct_value, cutoff_cycles, threshold_rfu):
    """
    Interprets a PCR result based on Ct value, cutoff, and threshold.

    Args:
        ct_value (float): The cycle threshold value.
        cutoff_cycles (int): The assay cutoff in cycles.
        threshold_rfu (int): The threshold in Relative Fluorescence Units.
    """
    print(f"Analyzing PCR Result:")
    print(f"  - Threshold: {threshold_rfu:,} RFU")
    print(f"  - Assay Cut-off: {cutoff_cycles} cycles")
    print(f"  - Sample Ct Value: {ct_value}")
    print("-" * 20)

    if ct_value < cutoff_cycles:
        interpretation = (
            f"The result is POSITIVE.\n"
            f"The sample's amplification curve crossed the threshold of {threshold_rfu:,} RFU at cycle {ct_value}, "
            f"which is before the assay cut-off of {cutoff_cycles} cycles.\n"
            f"A Ct value of {ct_value} suggests a relatively low concentration of the target nucleic acid."
        )
    else:
        interpretation = (
            f"The result is NEGATIVE.\n"
            f"The sample's amplification curve did not cross the threshold of {threshold_rfu:,} RFU "
            f"before the assay cut-off of {cutoff_cycles} cycles."
        )

    print(interpretation)

# Given values from the problem
threshold = 125000
cutoff = 40
ct = 30.8

interpret_pcr_result(ct, cutoff, threshold)