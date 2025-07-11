def interpret_pcr_result(ct_value, cutoff_cycles, threshold_rfu):
    """
    Interprets a PCR result based on the Ct value and assay cut-off.
    
    Args:
        ct_value (float): The cycle at which fluorescence crossed the threshold.
        cutoff_cycles (int): The maximum number of cycles for the assay.
        threshold_rfu (int): The fluorescence threshold in RFU.
    """
    print("--- PCR Assay Parameters ---")
    print(f"Fluorescence Threshold: {threshold_rfu} RFU")
    print(f"Assay Cut-off: {cutoff_cycles} cycles")
    print(f"Measured Ct Value: {ct_value}")
    print("\n--- Interpretation Logic ---")
    print(f"The result is considered positive if the Ct value ({ct_value}) is less than the assay cut-off ({cutoff_cycles}).")
    
    # Check if the Ct value is below the cut-off
    if ct_value < cutoff_cycles:
        result = "POSITIVE"
        explanation = (
            f"The amplification curve crossed the threshold at cycle {ct_value}, "
            f"which is before the cut-off of {cutoff_cycles} cycles."
        )
    else:
        result = "NEGATIVE"
        explanation = (
            f"The amplification curve did not cross the threshold before the "
            f"cut-off of {cutoff_cycles} cycles."
        )
        
    print("\n--- Final Result ---")
    print(f"The sample is {result}.")
    print(f"Reasoning: {explanation}")

# Given values from the problem
threshold = 125000
cutoff = 40
ct = 30.8

# Run the interpretation
interpret_pcr_result(ct, cutoff, threshold)