def interpret_pcr_result(ct_value, assay_cutoff, threshold_rfu):
    """
    Interprets a PCR result based on the Ct value and assay cut-off.

    Args:
        ct_value (float): The cycle at which fluorescence crossed the threshold.
        assay_cutoff (int): The maximum number of cycles for the assay.
        threshold_rfu (int): The fluorescence threshold in RFU.
    """
    print("--- PCR Assay Interpretation ---")
    print(f"Fluorescence Threshold: {threshold_rfu} RFU")
    print(f"Assay Cycle Cut-off: {assay_cutoff} cycles")
    print(f"Measured Ct Value: {ct_value}")
    print("--------------------------------\n")
    
    print("Logic: Is the Ct value ({}) less than the assay cut-off ({})?".format(ct_value, assay_cutoff))

    if ct_value < assay_cutoff:
        result = "Positive"
        explanation = (
            f"Result: {result}\n"
            f"The result is interpreted as Positive because the amplification curve crossed the threshold\n"
            f"at cycle {ct_value}, which is before the assay cut-off of {assay_cutoff} cycles."
        )
    else:
        result = "Negative"
        explanation = (
            f"Result: {result}\n"
            f"The result is interpreted as Negative because the amplification curve did not cross the threshold\n"
            f"before the assay cut-off of {assay_cutoff} cycles."
        )
        
    print(explanation)
    return result

# Given parameters from the problem
pcr_threshold_rfu = 125000
pcr_assay_cutoff = 40
pcr_ct_value = 30.8

# Run the interpretation
final_interpretation = interpret_pcr_result(
    ct_value=pcr_ct_value,
    assay_cutoff=pcr_assay_cutoff,
    threshold_rfu=pcr_threshold_rfu
)