def interpret_pcr_result(ct_value, assay_cutoff):
    """
    Interprets a PCR result based on the Ct value and assay cut-off.

    Args:
        ct_value (float): The cycle threshold value.
        assay_cutoff (int): The maximum number of cycles for the assay.
    """
    print("--- PCR Interpretation ---")
    print(f"Assay Parameters:")
    print(f"  - Ct Value: {ct_value}")
    print(f"  - Assay Cut-off: {assay_cutoff} cycles\n")

    print("Logic:")
    print(f"The sample is 'Positive' if the Ct value ({ct_value}) is less than the assay cut-off ({assay_cutoff}).")
    print("-" * 28)

    if ct_value < assay_cutoff:
        result = "Positive"
        conclusion = f"Conclusion: The sample is {result} because the Ct value of {ct_value} is less than the assay cut-off of {assay_cutoff}."
    else:
        result = "Negative"
        conclusion = f"Conclusion: The sample is {result} because the Ct value of {ct_value} is not less than the assay cut-off of {assay_cutoff}."

    print(conclusion)
    return result

# Given values from the problem
given_ct_value = 30.8
given_assay_cutoff = 40

# Run the interpretation
interpret_pcr_result(given_ct_value, given_assay_cutoff)
