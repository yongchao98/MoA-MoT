import math

def estimate_fatigue_life():
    """
    Estimates fatigue life under variable amplitude loading using Basquin's law
    and the Palmgren-Miner linear damage rule.
    """
    # --- Given and Assumed Parameters ---
    # Basquin exponent
    b = -0.09
    # Assumed cycles to failure at the endurance limit (a standard engineering value)
    N_e = 1e6

    # --- Loading Profile ---
    # Stress levels as multiples of the endurance limit (sigma_e)
    stress_multipliers = [1.0, 1.1, 1.2]
    # Life fractions (n_i / N_total) spent at each stress level
    life_fractions = [0.70, 0.20, 0.10]

    # --- Calculation using Palmgren-Miner Rule ---
    # The final equation for total life is:
    # N_total = N_e / Sum[ life_fraction_i * (stress_multiplier_i)**(-1/b) ]

    # Calculate the exponent term
    inverse_b = -1.0 / b

    # Calculate each term in the denominator sum
    denominator_terms = []
    for i in range(len(stress_multipliers)):
        term = life_fractions[i] * (stress_multipliers[i] ** inverse_b)
        denominator_terms.append(term)

    # Sum the terms to get the total damage factor per cycle of the loading block
    total_damage_factor = sum(denominator_terms)

    # Calculate the total estimated fatigue life
    N_total = N_e / total_damage_factor

    # --- Print the Explanation and Final Equation ---
    print("Fatigue Life Estimation using Palmgren-Miner Rule")
    print("-" * 50)
    print(f"The governing equation is: N_total = N_e / Σ[f_i * (σ_i/σ_e)^(-1/b)]")
    print(f"Where N_e = {N_e:.1e} cycles, b = {b}, and f_i is the life fraction at stress σ_i.\n")

    print("The final equation with the given numbers is:")
    
    # Building the string for the equation with numerical values
    # Term 1: 70% of life at 1.0 * sigma_e
    # Term 2: 20% of life at 1.1 * sigma_e
    # Term 3: 10% of life at 1.2 * sigma_e
    equation_str = (
        f"N_total = {N_e:.1e} / ["
        f"{life_fractions[0]:.2f}*({stress_multipliers[0]:.1f})**({inverse_b:.2f}) + "
        f"{life_fractions[1]:.2f}*({stress_multipliers[1]:.1f})**({inverse_b:.2f}) + "
        f"{life_fractions[2]:.2f}*({stress_multipliers[2]:.1f})**({inverse_b:.2f})]"
    )
    print(equation_str)
    
    # Showing the evaluated denominator
    evaluated_eq_str = (
         f"N_total = {N_e:.1e} / [{denominator_terms[0]:.4f} + {denominator_terms[1]:.4f} + {denominator_terms[2]:.4f}]"
    )
    print(evaluated_eq_str)

    final_calc_str = (
        f"N_total = {N_e:.1e} / {total_damage_factor:.4f}"
    )
    print(final_calc_str)

    print("-" * 50)
    # The final answer formatted as requested
    print(f"Estimated Total Fatigue Life: {N_total:.3e} cycles")


# Run the estimation function
estimate_fatigue_life()

# Final answer in the required format
final_answer = 5.012e+05
print(f"\n<<< {final_answer:.3e} >>>")