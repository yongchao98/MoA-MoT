def compare_toxicity():
    """
    Compares the acute toxicity of TMT-Cl and TBT-Cl using their LD50 values.
    """
    # LD50 data for oral administration in rats (values can vary slightly between studies)
    tmt_cl_ld50 = 12.6  # mg/kg
    tbt_cl_ld50 = 129    # mg/kg

    print("--- Comparing Acute Toxicity of Organotin Compounds ---")
    print("The danger of a substance is often quantified by its LD50 (Lethal Dose, 50%).")
    print("A lower LD50 value indicates a higher level of toxicity.\n")

    print("Compound 1: Trimethyltin chloride (TMT-Cl)")
    print(f"Oral LD50 in rats: {tmt_cl_ld50} mg/kg")

    print("\nCompound 2: Tributyltin chloride (TBT-Cl)")
    print(f"Oral LD50 in rats: {tbt_cl_ld50} mg/kg")

    print("\n--- Analysis ---")
    # We are comparing the numbers 12.6 and 129
    print(f"The LD50 of TMT-Cl ({tmt_cl_ld50}) is significantly lower than the LD50 of TBT-Cl ({tbt_cl_ld50}).")

    if tmt_cl_ld50 < tbt_cl_ld50:
        factor = tbt_cl_ld50 / tmt_cl_ld50
        print(f"This means TMT-Cl is approximately {factor:.1f} times more acutely toxic than TBT-Cl.")
    
    print("\nThis large, measurable difference in lethality (LD50) is the most critical factor proving that TMT-Cl is more dangerous for humans.")

# Run the comparison
compare_toxicity()