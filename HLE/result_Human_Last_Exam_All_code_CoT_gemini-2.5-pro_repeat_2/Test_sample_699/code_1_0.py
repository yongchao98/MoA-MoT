import math

def estimate_fatigue_life():
    """
    Estimates the fatigue life of a specimen using the Basquin relationship
    and the Palmgren-Miner linear damage rule.
    """
    # Basquin exponent
    b = -0.09

    # Stress levels are given as multiples of the fatigue endurance limit (sigma_e)
    stress_ratios = [1.0, 1.1, 1.2]

    # Percentage of life spent at each stress level, converted to fractions
    cycle_fractions = [0.70, 0.20, 0.10]

    print("Step 1: Calculate the number of cycles to failure (N_f) for each stress level.")
    print("-----------------------------------------------------------------------------")
    print(f"Using Basquin's relation: N_f = (sigma_a / C)^(1/b)")
    print(f"With C = (2 * sigma_e) / (0.5^b), this simplifies to: N_f = ((ratio * (0.5^b)) / 2)^(1/b)\n")

    # Calculate N_f for each stress level
    N_f_values = []
    for ratio in stress_ratios:
        # N_f = ( (ratio * sigma_e) / C )^(1/b)
        # N_f = ( (ratio * sigma_e) / ( (2*sigma_e) / (0.5**b) ) )**(1/b)
        # N_f = ( (ratio * (0.5**b)) / 2 )**(1/b)
        n_f = ((ratio * (0.5**b)) / 2)**(1/b)
        N_f_values.append(n_f)

    (N_f1, N_f2, N_f3) = N_f_values
    print(f"Fatigue life at {stress_ratios[0]}*sigma_e (N_f1): {N_f1:.3f} cycles")
    print(f"Fatigue life at {stress_ratios[1]}*sigma_e (N_f2): {N_f2:.3f} cycles")
    print(f"Fatigue life at {stress_ratios[2]}*sigma_e (N_f3): {N_f3:.3f} cycles\n")

    print("Step 2: Apply the Palmgren-Miner linear damage rule.")
    print("---------------------------------------------------------")
    print("Total Life (N_life) = 1 / sum(f_i / N_fi)\n")

    # Calculate the total damage sum in the denominator
    damage_sum = (cycle_fractions[0] / N_f1) + \
                 (cycle_fractions[1] / N_f2) + \
                 (cycle_fractions[2] / N_f3)

    # Calculate the total fatigue life
    total_life = 1 / damage_sum

    print("Final Equation:")
    final_equation = f"N_life = 1 / ({cycle_fractions[0]:.2f}/{N_f1:.3f} + {cycle_fractions[1]:.2f}/{N_f2:.3f} + {cycle_fractions[2]:.2f}/{N_f3:.3f})"
    print(final_equation)
    print(f"\nTotal Estimated Fatigue Life = {total_life:.3f} cycles")
    
    # Format the final answer in scientific notation
    formatted_answer = f"{total_life:.3e}"
    print(f"Final Answer (in Aeb format): {formatted_answer}")
    
    # The final answer tag as requested
    print(f"\n<<<{formatted_answer}>>>")


if __name__ == "__main__":
    estimate_fatigue_life()