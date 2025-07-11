import sys

# Disable the creation of __pycache__
sys.dont_write_bytecode = True

def interpret_pcr_result(ct_value, cutoff_cycles, threshold_rfu):
    """
    Interprets a PCR result based on the Ct value and assay cut-off.

    Args:
        ct_value (float): The cycle threshold value.
        cutoff_cycles (int): The assay cut-off in cycles.
        threshold_rfu (int): The fluorescence threshold in RFU.
    """
    print("--- PCR Assay Interpretation ---")
    print(f"Fluorescence Threshold: {threshold_rfu:,} RFU")
    print(f"Assay Cut-off: {cutoff_cycles} cycles")
    print(f"Measured Ct Value: {ct_value} cycles")
    print("-" * 35)
    
    # The main logic for interpretation
    is_positive = ct_value is not None and ct_value < cutoff_cycles
    
    print("Rule: A result is 'Positive' if the Ct value is less than the assay cut-off.")
    print(f"Evaluation: Is {ct_value} < {cutoff_cycles}?")
    print("-" * 35)

    if is_positive:
        print("Conclusion: The result is POSITIVE.")
        print("\nReasoning: The amplification curve crossed the fluorescence threshold")
        print(f"at cycle {ct_value}, which is before the cut-off of {cutoff_cycles} cycles.")
        print("This indicates that the target genetic material was detected.")
    else:
        print("Conclusion: The result is NEGATIVE.")
        print("\nReasoning: The amplification signal did not cross the threshold")
        print(f"before the cut-off of {cutoff_cycles} cycles.")
        print("This indicates that the target genetic material was not detected.")

if __name__ == '__main__':
    # Parameters from the provided information
    pcr_threshold_rfu = 125000
    pcr_assay_cutoff = 40
    pcr_ct_value = 30.8
    
    interpret_pcr_result(pcr_ct_value, pcr_assay_cutoff, pcr_threshold_rfu)
