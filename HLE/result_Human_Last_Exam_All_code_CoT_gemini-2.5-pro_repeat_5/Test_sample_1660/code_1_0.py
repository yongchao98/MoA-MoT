def compare_toxicity():
    """
    Compares the acute toxicity of Trimethyltin chloride (TMT-Cl) and
    Tributyltin chloride (TBT-Cl) based on their LD50 values.
    """
    # LD50 (oral, rat) values are a standard measure of acute toxicity.
    # A lower LD50 value means higher toxicity.
    ld50_tmt_cl_mg_per_kg = 12.6
    ld50_tbt_cl_mg_per_kg = 129.0

    print("Comparing the acute toxicity of two compounds based on LD50 values.")
    print(f"Trimethyltin chloride (TMT-Cl) LD50: {ld50_tmt_cl_mg_per_kg} mg/kg")
    print(f"Tributyltin chloride (TBT-Cl) LD50: {ld50_tbt_cl_mg_per_kg} mg/kg")
    print("-" * 30)

    # The final 'equation' is the comparison of these two numbers.
    print(f"Final Comparison Equation: {ld50_tmt_cl_mg_per_kg} < {ld50_tbt_cl_mg_per_kg}")

    if ld50_tmt_cl_mg_per_kg < ld50_tbt_cl_mg_per_kg:
        print("\nConclusion: TMT-Cl has a significantly lower LD50, meaning it is substantially more toxic and therefore more dangerous than TBT-Cl.")
    else:
        print("\nConclusion: TBT-Cl has a lower or equal LD50, meaning it is more or equally dangerous than TMT-Cl.")

compare_toxicity()