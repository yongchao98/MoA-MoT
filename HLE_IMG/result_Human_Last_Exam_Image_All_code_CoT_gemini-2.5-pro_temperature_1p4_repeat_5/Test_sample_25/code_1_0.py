import sys

def interpret_pcr_result(threshold, cutoff, ct_value):
    """
    Interprets a PCR assay result based on threshold, cutoff, and Ct value.
    Prints a step-by-step explanation and the final conclusion.
    """
    
    # Redirect print to a string to capture the final answer without printing it to the console with the explanation
    original_stdout = sys.stdout 
    
    print("--- PCR Result Interpretation ---")
    print(f"1. Assay Parameters:")
    print(f"   - Fluorescence Threshold: {int(threshold):,} RFU")
    print(f"   - Assay Cut-off: {int(cutoff)} cycles")
    print(f"   - Measured Ct Value: {ct_value} cycles")
    print("\n2. Interpretation Logic:")
    print("   The result is considered 'Positive' if the amplification signal crosses the threshold")
    print("   before the assay cut-off is reached. This is checked by comparing the Ct value to the cut-off.")
    
    # This represents the final logical check or 'equation'
    print("\n3. Evaluation:")
    print(f"   Is the Ct value ({ct_value}) less than the cut-off ({int(cutoff)})?")

    if ct_value is not None and ct_value < cutoff:
        result = "Positive"
        conclusion = f"   Yes, {ct_value} is less than {int(cutoff)}."
    else:
        result = "Negative"
        conclusion = f"   No, {ct_value} is not less than {int(cutoff)}."
    
    print(conclusion)
    
    print("\n--- Conclusion ---")
    print(f"The result is interpreted as: {result}")
    print("This indicates that the target genetic material was successfully detected and amplified.")
    
    # Final answer for parsing
    sys.stdout = original_stdout # Restore stdout
    return f"<<<{result}>>>"


# Given parameters from the problem description
threshold_rfu = 125000
assay_cutoff_cycles = 40
ct_value = 30.8

# Run the interpretation and capture the final answer
final_answer = interpret_pcr_result(threshold_rfu, assay_cutoff_cycles, ct_value)
print(final_answer)
