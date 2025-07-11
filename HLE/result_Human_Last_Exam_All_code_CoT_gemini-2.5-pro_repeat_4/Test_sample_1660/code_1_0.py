def compare_toxicity():
    """
    Compares the acute toxicity of Trimethyltin chloride (TMT-Cl) and
    Tributyltin chloride (TBT-Cl) based on their LD50 values.
    """
    # LD50 values (oral, rat) in mg/kg. Lower value means higher toxicity.
    ld50_tmt_cl = 13  # Approximate LD50 for Trimethyltin chloride
    ld50_tbt_cl = 130 # Approximate LD50 for Tributyltin chloride

    print("Comparing the acute toxicity of two compounds based on LD50 values.")
    print(f"LD50 is the dose (in mg/kg) that is lethal to 50% of a test population.")
    print("A lower LD50 value indicates higher toxicity.\n")

    print(f"Trimethyltin chloride (TMT-Cl) LD50: {ld50_tmt_cl} mg/kg")
    print(f"Tributyltin chloride (TBT-Cl) LD50: {ld50_tbt_cl} mg/kg\n")

    # The ratio shows how many times more toxic TMT-Cl is than TBT-Cl in terms of acute lethality.
    toxicity_ratio = ld50_tbt_cl / ld50_tmt_cl

    if ld50_tmt_cl < ld50_tbt_cl:
        print(f"Conclusion: TMT-Cl has a significantly lower LD50 value than TBT-Cl.")
        print(f"Based on these values, TMT-Cl is approximately {toxicity_ratio:.1f} times more acutely toxic.")
        print("This is the most direct and important factor making TMT-Cl more dangerous than TBT-Cl.")
    else:
        print("Conclusion: TBT-Cl has a lower or equal LD50 value, making it more or equally toxic.")

compare_toxicity()