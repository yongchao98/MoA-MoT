import math

def estimate_fatigue_life():
    """
    Estimates fatigue life based on the Basquin relationship and Palmgren-Miner rule.
    """
    # Given parameters from the problem description
    b = -0.09  # Basquin exponent

    # Life fractions (n_i / N_L) and their corresponding stress ratios (r_i)
    fractions = [0.70, 0.20, 0.10]
    stress_ratios = [1.0, 1.1, 1.2]

    # The endurance limit life, N_fe, is not given.
    # A standard engineering convention is to assume N_fe = 1,000,000 cycles.
    N_fe = 1e6

    print("Step 1: Define the problem parameters and assumptions.")
    print(f"Basquin exponent (b): {b}")
    print(f"Assumed life at endurance limit (N_fe): {int(N_fe)} cycles")
    print("-" * 30)

    # The governing equation derived from the Palmgren-Miner rule is:
    # N_L = N_fe / D_factor
    # where D_factor = sum(fraction_i * stress_ratio_i ^ (-1/b))

    # Calculate the exponent used in the damage factor calculation
    exponent = -1 / b

    print("Step 2: Calculate the combined damage factor.")
    print("Damage Factor = f1*r1^(-1/b) + f2*r2^(-1/b) + f3*r3^(-1/b)")
    
    term1 = fractions[0] * (stress_ratios[0] ** exponent)
    term2 = fractions[1] * (stress_ratios[1] ** exponent)
    term3 = fractions[2] * (stress_ratios[2] ** exponent)
    
    damage_factor = term1 + term2 + term3

    # Output the components of the damage factor equation
    print(f"Term 1 (70% life at 1.0*sigma_e): {fractions[0]:.2f} * {stress_ratios[0]:.1f}^({exponent:.3f}) = {term1:.4f}")
    print(f"Term 2 (20% life at 1.1*sigma_e): {fractions[1]:.2f} * {stress_ratios[1]:.1f}^({exponent:.3f}) = {term2:.4f}")
    print(f"Term 3 (10% life at 1.2*sigma_e): {fractions[2]:.2f} * {stress_ratios[2]:.1f}^({exponent:.3f}) = {term3:.4f}")
    print(f"Total Damage Factor = {term1:.4f} + {term2:.4f} + {term3:.4f} = {damage_factor:.4f}")
    print("-" * 30)

    # Step 3: Calculate the final estimated fatigue life (N_L)
    N_L = N_fe / damage_factor
    
    print("Step 3: Calculate the total estimated fatigue life (N_L).")
    print(f"N_L = N_fe / Damage Factor")
    print(f"N_L = {int(N_fe)} / {damage_factor:.4f}")
    print(f"N_L = {N_L:.3f} cycles")
    print("-" * 30)
    
    # Format the final answer into the required Aeb format
    final_answer_str = f"{N_L:.3e}".replace('e+0', 'e')
    
    print("Final answer in the required format (Aeb with 3 decimal places):")
    print(f"<<<{final_answer_str}>>>")

# Execute the function
estimate_fatigue_life()