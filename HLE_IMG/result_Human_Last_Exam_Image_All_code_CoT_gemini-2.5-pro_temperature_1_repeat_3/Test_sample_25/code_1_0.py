def interpret_pcr_result(ct_value, cutoff_cycles, threshold_rfu):
    """
    Interprets the result of a PCR assay based on Ct value and cut-off.
    """
    print("--- PCR Result Interpretation ---")
    print(f"Threshold (RFU): {threshold_rfu}")
    print(f"Assay Cut-off (Cycles): {cutoff_cycles}")
    print(f"Measured Ct Value: {ct_value}")
    print("-" * 33)

    # The primary rule for a positive result in qPCR
    if 0 < ct_value < cutoff_cycles:
        print("\nConclusion: The result is POSITIVE.")
        print(f"The interpretation is based on the following comparison:")
        print(f"Ct Value ({ct_value}) < Assay Cut-off ({cutoff_cycles})")
        print("\nThis means the fluorescence signal crossed the threshold before the final cycle, indicating that the target nucleic acid was detected.")
    else:
        print("\nConclusion: The result is NEGATIVE.")
        print("The fluorescence signal did not cross the threshold before the assay cut-off cycle.")

# Given parameters
threshold = 125000
assay_cutoff = 40
ct = 30.8

# Run the interpretation
interpret_pcr_result(ct, assay_cutoff, threshold)