import math

def final_calculation():
    """
    This function performs the derivation and calculation as per the problem's constraints.
    """
    
    # Step 1: Define true values for final error calculation
    rho_true = 0.9
    r_true = 0.5
    pi_true = math.pi
    mass_true = rho_true * (4/3) * pi_true * (r_true**3)
    
    print("Derivation of the mass of the rock using Titan 5-bit fractional arithmetic.\n")
    print(f"The exact formula is: mass = density * (4/3) * pi * radius^3")
    print(f"The true mass is approximately {mass_true:.7f} kg.\n")

    # Step 2: Define inputs as 5-bit fractions
    # Numerators and denominators must be <= 31
    rho_num, rho_den = 9, 10
    r_cubed_num, r_cubed_den = 1, 8
    four_thirds_num, four_thirds_den = 4, 3
    
    # Crucial choice for pi to make the calculation possible without overflow
    pi_approx_num, pi_approx_den = 28, 9

    print("Step 1: Represent all values as 5-bit fractions.")
    print(f"density (rho) = {rho_num}/{rho_den}")
    print(f"4/3 = {four_thirds_num}/{four_thirds_den}")
    print(f"radius^3 (r^3) = (1/2)^3 = {r_cubed_num}/{r_cubed_den}")
    print(f"pi is approximated as: pi ≈ {pi_approx_num}/{pi_approx_den}\n")
    
    # Step 3: Write the full expression
    print("Step 2: The full expression for the mass is:")
    print(f"mass = ({rho_num}/{rho_den}) * ({four_thirds_num}/{four_thirds_den}) * ({pi_approx_num}/{pi_approx_den}) * ({r_cubed_num}/{r_cubed_den})\n")

    # Step 4: Perform the calculation step-by-step, grouping to avoid overflow
    print("Step 3: Calculate step-by-step, simplifying to keep numbers within 5-bit limits (<= 31).")
    print("We group terms strategically: mass = [ (9/10) * (28/9) ] * [ (4/3) * (1/8) ]\n")
    
    # First group: (9/10) * (28/9)
    # Simplify by cancelling 9s: 28/10, then reduce to 14/5
    group1_num, group1_den = 14, 5
    print(f"Calculating the first group: (9/10) * (28/9)")
    print(f"   The 9 in the numerator and denominator cancel out, leaving 28/10.")
    print(f"   28/10 is simplified to {group1_num}/{group1_den}. Both numbers are valid (<= 31).\n")
    
    # Second group: (4/3) * (1/8)
    # Simplify 4/8 to 1/2: 1/(3*2) = 1/6
    group2_num, group2_den = 1, 6
    print(f"Calculating the second group: (4/3) * (1/8)")
    print(f"   This becomes 4/24, which simplifies to {group2_num}/{group2_den}. Both numbers are valid (<= 31).\n")
    
    # Final multiplication: (14/5) * (1/6)
    # Becomes 14/30, which simplifies to 7/15
    final_num, final_den = 7, 15
    print(f"Final step, multiply the results of the two groups: ({group1_num}/{group1_den}) * ({group2_num}/{group2_den})")
    print(f"   This becomes 14/30, which simplifies to the final result: {final_num}/{final_den}.\n")
    
    # Step 5: Calculate final error
    mass_calc = final_num / final_den
    abs_error = abs(mass_true - mass_calc)
    
    print("Final Result:")
    print(f"The calculated mass is {final_num}/{final_den}, which is approximately {mass_calc:.7f} kg.")
    print(f"The absolute error is |{mass_true:.7f} - {mass_calc:.7f}| = {abs_error:.7f}")
    
    # The question asks for the error rounded to 0.001 (3 decimal places)
    final_error_rounded = round(abs_error, 3)
    
    print(f"The smallest absolute error, rounded to 0.001, is: {final_error_rounded}")
    print(f"<<<{final_error_rounded}>>>")

final_calculation()