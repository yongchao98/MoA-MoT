import math

def estimate_fatigue_life():
    """
    Estimates fatigue life using the Basquin relationship and Palmgren-Miner rule.
    """
    # Step 1: Define constants and assumptions
    b = -0.09  # Basquin exponent
    # Assume endurance limit corresponds to failure at N_fe cycles
    N_fe = 1e7

    # Stress levels are given as fractions of life at certain stress ratios
    # c_i: fraction of life at a given stress level
    # sr_i: stress ratio (sigma_i / sigma_e)
    c = [0.7, 0.2, 0.1]
    sr = [1.0, 1.1, 1.2]

    # Step 2: Calculate the denominator of the Palmgren-Miner equation
    # The total life N_total = N_fe / [c1*(sr1)^(-1/b) + c2*(sr2)^(-1/b) + c3*(sr3)^(-1/b)]
    # Let's calculate the value of each term in the denominator.
    # The exponent is -1/b
    exponent = -1.0 / b

    term1 = c[0] * math.pow(sr[0], exponent)
    term2 = c[1] * math.pow(sr[1], exponent)
    term3 = c[2] * math.pow(sr[2], exponent)
    
    denominator = term1 + term2 + term3

    # Step 3: Calculate the total fatigue life
    N_total = N_fe / denominator

    # Step 4: Print the full equation with intermediate values and the final result
    print("Fatigue Life Estimation using Palmgren-Miner Rule:")
    print("-" * 55)
    print(f"Basquin Exponent (b): {b}")
    print(f"Assumed Cycles at Endurance Limit (N_fe): {N_fe:.3e}")
    print(f"Exponent for calculation (-1/b): {-1/b:.3f}")
    print("\nPalmgren-Miner Equation for Total Life (N_total):")
    print("N_total = N_fe / [c1*(sr1)^(-1/b) + c2*(sr2)^(-1/b) + c3*(sr3)^(-1/b)]\n")

    print("Substituting the values:")
    # Print the equation with all the numbers
    print(f"N_total = {N_fe:.3e} / [{c[0]}*({sr[0]})^{exponent:.3f} + {c[1]}*({sr[1]})^{exponent:.3f} + {c[2]}*({sr[2]})^{exponent:.3f}]")
    print(f"N_total = {N_fe:.3e} / [{term1:.3f} + {term2:.3f} + {term3:.3f}]")
    print(f"N_total = {N_fe:.3e} / {denominator:.3f}")

    print("\nFinal Estimated Fatigue Life:")
    # The output format Aeb is scientific notation, represented by 'e' in python f-strings.
    print(f"{N_total:.3e} cycles")

    # Final answer in the requested format
    final_answer = f"{N_total:.3e}"
    print(f"\n<<<{final_answer}>>>")

if __name__ == '__main__':
    estimate_fatigue_life()