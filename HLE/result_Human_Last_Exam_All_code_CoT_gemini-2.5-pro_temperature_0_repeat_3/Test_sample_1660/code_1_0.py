def compare_toxicity():
    """
    Compares the toxicity of TMT-Cl and TBT-Cl based on their LD50 values.
    LD50 (Lethal Dose, 50%) is the amount of a substance that kills 50% of a test population.
    A lower LD50 value indicates higher toxicity.
    """
    # LD50 values (oral, rat) in mg/kg.
    ld50_tmt_cl = 12.6  # Trimethyltin chloride
    ld50_tbt_cl = 175   # Tributyltin chloride (an average value from the typical range)

    print("Toxicity Comparison:")
    print(f"Trimethyltin chloride (TMT-Cl) LD50: {ld50_tmt_cl} mg/kg")
    print(f"Tributyltin chloride (TBT-Cl) LD50: {ld50_tbt_cl} mg/kg")
    print("-" * 40)

    if ld50_tmt_cl < ld50_tbt_cl:
        # Calculate how many times more toxic TMT-Cl is based on LD50 values
        toxicity_factor = ld50_tbt_cl / ld50_tmt_cl
        print(f"The LD50 of TMT-Cl ({ld50_tmt_cl}) is significantly lower than that of TBT-Cl ({ld50_tbt_cl}).")
        print(f"This means TMT-Cl is approximately {toxicity_factor:.1f} times more acutely toxic.")
        print("\nConclusion: The most important factor making TMT-Cl more dangerous is its significantly lower LD50 value, which is a direct measure of its high acute toxicity.")
    else:
        # This case is not expected based on known toxicological data
        print("Based on the data, TBT-Cl appears more or equally toxic to TMT-Cl.")

compare_toxicity()