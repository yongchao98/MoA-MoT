def compare_toxicity():
    """
    Compares the toxicity of Tributyltin chloride (TBT-Cl) and
    Trimethyltin chloride (TMT-Cl) based on their LD50 values.
    """

    # LD50 values (oral, rat) in mg/kg.
    # LD50 (Lethal Dose, 50%) is the dose required to kill half the members of a tested population.
    # A lower LD50 value indicates higher toxicity.
    ld50_tmt_cl = 12.6  # Trimethyltin chloride
    ld50_tbt_cl = 129.0 # Tributyltin chloride

    print("Comparing the acute toxicity of two compounds:")
    print("-" * 50)
    print(f"Trimethyltin chloride (TMT-Cl) LD50: {ld50_tmt_cl} mg/kg")
    print(f"Tributyltin chloride (TBT-Cl) LD50:  {ld50_tbt_cl} mg/kg")
    print("-" * 50)

    # Compare the values and draw a conclusion
    if ld50_tmt_cl < ld50_tbt_cl:
        ratio = ld50_tbt_cl / ld50_tmt_cl
        print(f"Conclusion: TMT-Cl has a significantly lower LD50 value than TBT-Cl.")
        print(f"This means TMT-Cl is approximately {ratio:.1f} times more acutely toxic.")
        print("\nTherefore, the most important factor making TMT-Cl more dangerous is its significantly lower LD50 value.")
    else:
        print("Conclusion: TBT-Cl has a lower or equal LD50, which contradicts the premise.")

# Run the comparison
compare_toxicity()