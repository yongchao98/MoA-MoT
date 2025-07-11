import sys

def interpret_pcr_result(ct_value, cutoff_cycles, threshold_rfu):
    """
    Interprets a PCR result based on Ct value and assay cutoff.

    Args:
        ct_value (float): The cycle threshold value.
        cutoff_cycles (int): The assay cutoff in cycles.
        threshold_rfu (int): The fluorescence threshold in RFU.
    """
    print("--- PCR Assay Interpretation ---")
    print(f"Analysis based on the following parameters:")
    print(f"  - Fluorescence Threshold: {threshold_rfu:,} RFU")
    print(f"  - Assay Cut-off: {cutoff_cycles} cycles")
    print(f"  - Measured Ct Value: {ct_value}")
    print("----------------------------------")

    # A result is positive if the Ct value is valid and less than the cutoff.
    if ct_value is not None and ct_value < cutoff_cycles:
        result = "Positive"
        reason = (f"The result is interpreted as POSITIVE because the Ct value of {ct_value} "
                  f"is less than the assay cut-off of {cutoff_cycles} cycles. "
                  f"This indicates that the amplification signal successfully crossed the {threshold_rfu:,} RFU threshold within the allowed time.")
    else:
        result = "Negative"
        reason = (f"The result is interpreted as NEGATIVE because the amplification signal "
                  f"did not cross the {threshold_rfu:,} RFU threshold before the assay cut-off of {cutoff_cycles} cycles.")

    print(reason)
    # The final answer format is handled outside the function
    # but we can return the concise result.
    return result


if __name__ == '__main__':
    # Given parameters
    threshold = 125000
    cutoff = 40
    ct = 30.8

    # Get the interpretation
    final_result = interpret_pcr_result(ct_value=ct, cutoff_cycles=cutoff, threshold_rfu=threshold)
    
    # Writing final answer to a file to be captured by the system, as per instructions for some platforms.
    # For direct output, this part can be modified.
    # We will print it directly to stdout in this case.
    # The '<<<...>>>' format is a marker for the final answer.
    sys.stdout.write(f"\n<<<{final_result}>>>\n")
