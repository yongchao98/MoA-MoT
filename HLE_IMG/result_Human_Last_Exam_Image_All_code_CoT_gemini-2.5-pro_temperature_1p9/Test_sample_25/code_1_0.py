def interpret_pcr_result():
    """
    Interprets the result of a PCR assay based on threshold, cut-off, and Ct value.
    """
    # Define the parameters for the PCR assay interpretation
    threshold_rfu = 125000
    assay_cut_off_cycles = 40
    ct_value = 30.8

    # Print the parameters used for interpretation
    print("--- PCR Assay Interpretation ---")
    print(f"Fluorescence Threshold: {threshold_rfu:,} RFU")
    print(f"Assay Cut-off: {assay_cut_off_cycles} cycles")
    print(f"Measured Ct Value: {ct_value}")
    print("---------------------------------\n")

    # Determine the result based on the comparison
    print("Evaluation:")
    if ct_value < assay_cut_off_cycles:
        result_text = "POSITIVE"
        explanation = (f"The Ct value ({ct_value}) is less than the assay cut-off ({assay_cut_off_cycles}), "
                       "indicating that the target nucleic acid was detected.")
    else:
        result_text = "NEGATIVE"
        explanation = (f"The Ct value ({ct_value}) is not less than the assay cut-off ({assay_cut_off_cycles}), "
                       "indicating that the target nucleic acid was not detected.")
    
    print(f"Result: {result_text}")
    print(f"Reason: {explanation}")

# Run the interpretation function
interpret_pcr_result()
