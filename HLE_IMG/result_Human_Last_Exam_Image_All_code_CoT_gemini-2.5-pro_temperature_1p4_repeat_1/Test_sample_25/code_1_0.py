def interpret_pcr_result(ct_value, cutoff_cycles, threshold_rfu):
    """
    Interprets a PCR result based on the Ct value and assay cut-off.
    """
    print("--- PCR Assay Interpretation ---")
    print(f"Threshold RFU: {threshold_rfu}")
    print(f"Assay Cut-off: {cutoff_cycles} cycles")
    print(f"Measured Ct Value: {ct_value} cycles")
    print("-" * 30)
    print("Rule: The result is 'Positive' if the Ct value is less than the assay cut-off.")
    print("-" * 30)

    if ct_value < cutoff_cycles:
        result = "Positive"
        explanation = (
            f"Conclusion: The sample is POSITIVE. \n"
            f"The amplification curve crossed the threshold at cycle {ct_value}, "
            f"which occurred before the assay cut-off of {cutoff_cycles} cycles."
        )
    else:
        result = "Negative"
        explanation = (
            f"Conclusion: The sample is NEGATIVE. \n"
            f"The amplification curve did not cross the threshold before the "
            f"assay cut-off of {cutoff_cycles} cycles."
        )
    
    print(explanation)
    return result

# Given parameters from the problem description
pcr_threshold_rfu = 125000
pcr_cutoff_cycles = 40
pcr_ct_value = 30.8

# Run the interpretation
final_result = interpret_pcr_result(pcr_ct_value, pcr_cutoff_cycles, pcr_threshold_rfu)
