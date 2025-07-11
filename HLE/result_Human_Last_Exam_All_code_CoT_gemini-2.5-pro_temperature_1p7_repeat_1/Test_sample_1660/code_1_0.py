def toxicity_comparison():
    """
    Compares the acute toxicity of Trimethyltin chloride (TMT-Cl)
    and Tributyltin chloride (TBT-Cl) using their LD50 values.
    """

    # Define chemical data: name and LD50 value (oral, mouse) in mg/kg.
    tmt_cl_name = "Trimethyltin chloride (TMT-Cl)"
    tmt_cl_ld50 = 3  # A lower LD50 value means higher toxicity.

    tbt_cl_name = "Tributyltin chloride (TBT-Cl)"
    tbt_cl_ld50 = 129

    print("--- Comparing Chemical Toxicity using LD50 ---")
    print("LD50 is the dose required to be lethal to 50% of a test population.")
    print("A lower LD50 value indicates a more dangerous substance.\n")

    print(f"Substance: {tmt_cl_name}")
    print(f"LD50 value (oral, mouse): {tmt_cl_ld50} mg/kg\n")

    print(f"Substance: {tbt_cl_name}")
    print(f"LD50 value (oral, mouse): {tbt_cl_ld50} mg/kg\n")

    # Determine which chemical is more toxic and by what factor.
    if tmt_cl_ld50 < tbt_cl_ld50:
        more_toxic_chemical = tmt_cl_name
        less_toxic_chemical = tbt_cl_name
        factor = tbt_cl_ld50 / tmt_cl_ld50
    else:
        more_toxic_chemical = tbt_cl_name
        less_toxic_chemical = tmt_cl_name
        factor = tmt_cl_ld50 / tbt_cl_ld50

    print("--- Conclusion ---")
    print(f"The LD50 of {tmt_cl_name} ({tmt_cl_ld50} mg/kg) is significantly lower than that of {tbt_cl_name} ({tbt_cl_ld50} mg/kg).")
    print(f"This means {more_toxic_chemical} is approximately {factor:.0f} times more acutely toxic.")
    print("\nThis substantial difference in LD50 values is the most important factor making TMT-Cl more dangerous for human beings than TBT-Cl.")

# Execute the function to display the comparison.
toxicity_comparison()