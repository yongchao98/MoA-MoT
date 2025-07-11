def analyze_tree_ring_data():
    """
    This function analyzes the factors influencing the 13C isotope ratio
    in tree rings from 1886-1990 and determines the most likely cause
    from a given set of options.
    """
    
    print("Analyzing the predominant factor for the declining 13C ratio in Chinese pine tree rings (1886-1990 AD).")
    print("-" * 70)

    # Step 1: Explain the primary global cause, which is not listed as an option.
    print("Step 1: The 'Suess Effect'")
    print("The single most dominant factor causing a global decline in the 13C/12C ratio (δ¹³C) in tree rings since the Industrial Revolution is the burning of fossil fuels.")
    print("Fossil fuels are ancient plant matter, naturally low in 13C. Burning them releases CO2 with a low 13C ratio into the atmosphere.")
    print("This changes the isotopic composition of the air that trees use for growth, and it is the true predominant factor. However, it is not an answer choice.")
    print("-" * 70)

    # Step 2: Evaluate the provided options.
    print("Step 2: Evaluating the given answer choices.")
    print("\n[B] Periods of drought: Incorrect.")
    print("   Drought causes water stress. To save water, trees partially close their leaf pores (stomata), which leads to an INCREASE in the 13C ratio, the opposite of the observed trend.")

    print("\n[A, C, D] Other factors: Less likely to be predominant.")
    print("   Tree maturity (A), starch reserves (C), and earlywood proportion (D) are secondary or localized factors and are not powerful enough to be the main driver of a consistent, 100-year regional trend.")

    print("\n[E] Changes in the SE Asia monsoon: Most plausible option.")
    print("   This is a large-scale climate system. A long-term trend towards a stronger, wetter monsoon reduces water stress. With ample water, trees can keep their stomata fully open, which allows them to be more selective and take up less of the heavy 13C isotope. This physiological response leads to a DECLINING 13C ratio, which matches the observation.")
    print("-" * 70)
    
    # Step 3: Present the conclusion and a conceptual equation.
    print("Step 3: Conclusion and Illustrative Equation.")
    print("Given the options, 'Changes in the SE Asia monsoon' is the most plausible large-scale environmental driver.")
    
    print("\nA conceptual equation shows how different factors contribute to the final 13C value. A stronger monsoon contributes a negative value (a decline).")
    
    effect_monsoon = -0.4
    effect_drought_years = 0.2
    effect_local = -0.1
    total_climate_effect = effect_monsoon + effect_drought_years + effect_local
    
    print("Equation: Δ_13C_Climate = Effect_Monsoon + Effect_Drought + Effect_Local")
    print(f"Example: Δ_13C_Climate = {effect_monsoon} + {effect_drought_years} + {effect_local} = {total_climate_effect:.1f} (in per mil ‰)")

if __name__ == '__main__':
    analyze_tree_ring_data()