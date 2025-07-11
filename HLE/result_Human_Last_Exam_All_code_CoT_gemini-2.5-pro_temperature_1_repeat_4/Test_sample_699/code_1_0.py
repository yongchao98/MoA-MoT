import math

def estimate_fatigue_life():
    """
    Estimates fatigue life using the Palmgren-Miner rule and Basquin relationship.
    """
    # 1. Define constants and loading profile
    b = -0.09  # Basquin exponent
    N_e = 1e6  # Assumed cycles at endurance limit (a standard value)
    
    # Life fractions (n_i / L) for each stress level
    life_fractions = [0.70, 0.20, 0.10]
    
    # Stress ratios (sigma_i / sigma_e) for each stress level
    stress_ratios = [1.0, 1.1, 1.2]

    # 2. Calculate the denominator of the life equation
    # Denominator = Sum[ life_fraction_i / (stress_ratio_i)^(1/b) ]
    # This comes from rearranging the Palmgren-Miner rule: L = N_e / Sum(...)
    try:
        term1 = life_fractions[0] / (stress_ratios[0]**(1/b))
        term2 = life_fractions[1] / (stress_ratios[1]**(1/b))
        term3 = life_fractions[2] / (stress_ratios[2]**(1/b))
        denominator = term1 + term2 + term3
    except ZeroDivisionError:
        print("Error: The Basquin exponent 'b' cannot be zero.")
        return

    # 3. Calculate the total fatigue life L
    total_life = N_e / denominator

    # 4. Print the final equation with all numbers as requested
    print("The fatigue life (L) is calculated using the following equation based on the Palmgren-Miner rule:")
    
    # Construct and print the equation string with unevaluated powers
    equation_str = (
        f"L = {N_e:g} / ("
        f"{life_fractions[0]:.3f} / {stress_ratios[0]:.1f}^(1/{b}) + "
        f"{life_fractions[1]:.3f} / {stress_ratios[1]:.1f}^(1/{b}) + "
        f"{life_fractions[2]:.3f} / {stress_ratios[2]:.1f}^(1/{b}))"
    )
    print(equation_str)

    # Print the equation with evaluated terms for clarity
    evaluated_terms_str = (
         f"L = {N_e:g} / ({term1:.3f} + {term2:.3f} + {term3:.3f})"
    )
    print("\n" + evaluated_terms_str)
    
    # Print the final calculation step
    final_calc_str = f"L = {N_e:g} / {denominator:.3f}"
    print(final_calc_str)
    
    # 5. Print the final answer in the required format
    result_str = f"\nEstimated Fatigue Life (L) = {total_life:.3e} cycles"
    print(result_str)

# Execute the function to get the result
estimate_fatigue_life()
<<<4.849e+05>>>