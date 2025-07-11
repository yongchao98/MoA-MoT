def compare_toxicity():
    """
    Compares the toxicity of Trimethyltin chloride (TMT-Cl) and
    Tributyltin chloride (TBT-Cl) based on their LD50 values.
    """
    # LD50 values (oral, rat) in mg/kg are a standard measure of acute toxicity.
    # Data is sourced from toxicological databases.
    ld50 = {
        'TMT-Cl': 13,
        'TBT-Cl': 129
    }

    tmt_ld50 = ld50['TMT-Cl']
    tbt_ld50 = ld50['TBT-Cl']

    print("--- Comparing Chemical Toxicity ---")
    print("A lower LD50 value indicates higher acute toxicity.")
    print(f"LD50 of Trimethyltin chloride (TMT-Cl): {tmt_ld50} mg/kg")
    print(f"LD50 of Tributyltin chloride (TBT-Cl): {tbt_ld50} mg/kg")
    print("-" * 35)

    if tmt_ld50 < tbt_ld50:
        toxicity_ratio = tbt_ld50 / tmt_ld50
        print(f"TMT-Cl has a significantly lower LD50 value than TBT-Cl.")
        print(f"Based on this data, TMT-Cl is approximately {toxicity_ratio:.1f} times more acutely toxic.")
        print("\nConclusion:")
        print("The most important factor making TMT-Cl more dangerous is its significantly lower LD50 value,")
        print("which is a direct measure of its high acute toxicity.")
    else:
        print("The provided LD50 values suggest TBT-Cl is more toxic, which contradicts established data.")

compare_toxicity()