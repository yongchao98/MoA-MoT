import math

def estimate_fatigue_life():
    """
    Estimates the fatigue life of a specimen using the Palmgren-Miner rule
    and the Basquin relationship for the given loading conditions.
    """
    # Given parameters
    b = -0.09  # Basquin exponent
    
    # Life fractions and corresponding stress multipliers
    fractions = [0.70, 0.20, 0.10]
    stress_multipliers = [1.0, 1.1, 1.2]
    
    # The fatigue life (L) is given by the Palmgren-Miner rule:
    # L * [ (f1/N1) + (f2/N2) + (f3/N3) ] = 1
    # Where Ni is the life at stress level i.
    # From Basquin's law, Ni = Ne * (stress_multiplier_i)^(1/b)
    # The equation for the life ratio L/Ne becomes:
    # L/Ne = 1 / [ f1*(r1)^(-1/b) + f2*(r2)^(-1/b) + f3*(r3)^(-1/b) ]
    
    # The exponent in the damage sum
    inv_b_exp = -1.0 / b

    # Calculate each term in the denominator of the life ratio equation
    term1 = fractions[0] * (stress_multipliers[0] ** inv_b_exp)
    term2 = fractions[1] * (stress_multipliers[1] ** inv_b_exp)
    term3 = fractions[2] * (stress_multipliers[2] ** inv_b_exp)
    
    # Calculate the total denominator (which represents the total damage for one block of cycles)
    damage_sum_denominator = term1 + term2 + term3
    
    # The final fatigue life is the reciprocal of the damage sum
    fatigue_life_ratio = 1.0 / damage_sum_denominator

    # --- Output Section ---
    print("The fatigue life ratio (L/N_e) is calculated using the Palmgren-Miner rule:")
    print(f"L/N_e = 1 / ( f1*r1**(-1/b) + f2*r2**(-1/b) + f3*r3**(-1/b) )")
    print("\nSubstituting the given values into the equation:")
    print(f"Basquin exponent, b = {b}")
    print(f"Inverse exponent, -1/b = {inv_b_exp:.3f}")
    
    print("\nFinal equation with calculated components:")
    # We explicitly calculate the powered terms to show each number in the final equation
    r2_pow_val = stress_multipliers[1]**inv_b_exp
    r3_pow_val = stress_multipliers[2]**inv_b_exp
    
    print(f"Life Ratio = 1 / ( {fractions[0]:.2f}*{stress_multipliers[0]:.1f}**{inv_b_exp:.3f} + "
          f"{fractions[1]:.2f}*{stress_multipliers[1]:.1f}**{inv_b_exp:.3f} + "
          f"{fractions[2]:.2f}*{stress_multipliers[2]:.1f}**{inv_b_exp:.3f} )")

    print(f"Life Ratio = 1 / ( {fractions[0]:.3f} + {fractions[1]:.3f}*{r2_pow_val:.3f} + {fractions[2]:.3f}*{r3_pow_val:.3f} )")
    print(f"Life Ratio = 1 / ( {term1:.3f} + {term2:.3f} + {term3:.3f} )")
    print(f"Life Ratio = 1 / ( {damage_sum_denominator:.3f} )")

    print(f"\nEstimated Fatigue Life Ratio = {fatigue_life_ratio:.3e}")

# Execute the function
estimate_fatigue_life()
<<<4.625e-01>>>