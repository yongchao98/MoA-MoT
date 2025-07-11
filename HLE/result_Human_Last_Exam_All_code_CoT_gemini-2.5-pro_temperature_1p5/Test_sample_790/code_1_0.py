import math

def titan_calculation():
    """
    This function simulates the calculation of the required force on the Titan
    superconducting computer, following its specific computational rules.
    """
    print("### Titan Computer: Force Calculation for Coconut Problem ###\n")

    # 1. Physics Formulation
    # The rock is under a constant thrust F at 45 degrees and gravity mg downwards.
    # The net acceleration components are:
    # ax = (F * cos(45)) / m
    # ay = (F * sin(45) - m*g) / m
    # The trajectory is y/x = ay/ax. For a 45-degree angle, sin(45)=cos(45).
    # y/x = 1 - (m*g) / (F*cos(45))
    # For target (x=20, y=10), y/x = 1/2.
    # 1/2 = 1 - (m*g*sqrt(2)) / F
    # This gives the formula for the force F:
    # F = 2 * m * g * sqrt(2)
    print("Step 1: Derived the physics equation for the force F.")
    print("   F = 2 * m * g * sqrt(2)\n")

    # 2. Mass Calculation
    # r = 0.5 cm = 1/200 m
    # rho = 0.9 kg/cm^3 = 0.9 * (100^3) kg/m^3 = 900,000 kg/m^3
    # V = (4/3) * pi * r^3
    # m = rho * V = 900,000 * (4/3) * pi * (1/200)^3
    # m = 900,000 * (4/3) * pi * (1 / 8,000,000)
    # m = (900,000 / 8,000,000) * (4/3) * pi
    # m = (9/80) * (4/3) * pi = (3/20) * pi
    print("Step 2: Calculated the mass 'm' of the rock.")
    print("   m = (density) * (volume) = (900,000 kg/m^3) * (4/3 * pi * (0.005 m)^3)")
    print("   m = (3/20) * pi kg\n")

    # 3. Final Equation for F
    # F = 2 * ((3/20) * pi) * g * sqrt(2)
    # F = (3/10) * pi * g * sqrt(2)
    print("Step 3: Substituted 'm' into the force equation.")
    print("   F = (3/10) * pi * g * sqrt(2)\n")

    # 4. Titan Computation
    # We must use 5-bit integer fractions (numerators/denominators <= 31).
    # We choose approximations that allow for simplification.
    pi_approx = "22/7"
    g_approx = "10/1"
    sqrt2_approx = "7/5"
    print("Step 4: Performing calculation using Titan's 5-bit fractional arithmetic.")
    print(f"   Selected approximations:")
    print(f"   - pi     ≈ {pi_approx}")
    print(f"   - g      ≈ {g_approx}  (True value is ~9.8)")
    print(f"   - sqrt(2) ≈ {sqrt2_approx} (True value is ~1.414)\n")

    print("   Calculation Narrative:")
    # F = (3/10) * (22/7) * (10/1) * (7/5)
    # Group terms to perform cancellation.
    print("   F = ( (3/10) * (10/1) ) * ( (22/7) * (7/5) )")

    # Step 4a: First multiplication
    # (3/10) * (10/1) = (3*10)/(10*1) -> simplifies to 3/1
    term1_num, term1_den = 3, 1
    print(f"   - Multiplying (3/10) by (10/1). Cross-cancellation of 10 gives: {term1_num}/{term1_den}")

    # Step 4b: Second multiplication
    # (22/7) * (7/5) = (22*7)/(7*5) -> simplifies to 22/5
    term2_num, term2_den = 22, 5
    print(f"   - Multiplying (22/7) by (7/5). Cross-cancellation of 7 gives: {term2_num}/{term2_den}")

    # Step 4c: Final multiplication
    # Now we must compute (3/1) * (22/5).
    # The direct result would be (3*22)/5 = 66/5.
    # However, 66 is larger than 31 and cannot be stored in a 5-bit register.
    # An approximation is required as per Titan's rules.
    print(f"   - Next operation: ({term1_num}/{term1_den}) * ({term2_num}/{term2_den})")
    print(f"   - Numerator {term1_num}*{term2_num} = 66, which exceeds the 31 limit. Overflow!")
    
    # We approximate one of the terms to allow for cancellation.
    # The value of 22/5 is 4.4. A close fraction that has 3 in the denominator is 13/3 (~4.33).
    approx_num, approx_den = 13, 3
    print(f"   - Applying constraint rule: Approximate {term2_num}/{term2_den} (value 4.4) to {approx_num}/{approx_den} (~4.33) to enable simplification.")

    # Perform final calculation with the approximated fraction
    # (3/1) * (13/3) = (3*13)/(1*3) -> simplifies to 13/1
    final_num, final_den = 13, 1
    print(f"   - Final multiplication: ({term1_num}/{term1_den}) * ({approx_num}/{approx_den}) gives {final_num}/{final_den} after cancellation.\n")

    # 5. Final Result and Error
    final_force_val = final_num / final_den
    print("Step 5: Final Result and Error Analysis.")
    print("   The complete calculation is:")
    print(f"   F = (3/10) * pi * g * sqrt(2)")
    print(f"     ≈ (3/10) * (10/1) * (22/7) * (7/5)")
    print(f"     = (3/1) * (22/5)")
    print(f"     ≈ (3/1) * (13/3)  <- Approximation to satisfy 5-bit constraint")
    print(f"     = {final_num}/{final_den}")

    true_F = (3.0/10.0) * math.pi * 9.8 * math.sqrt(2.0)
    error = abs(final_force_val - true_F)

    print(f"\n   The force calculated by Titan is {final_force_val} N.")
    print(f"   The true theoretical force is approximately {true_F:.4f} N.")
    print(f"   The absolute error 'e' is |{final_force_val} - {true_F:.4f}| = {error:.4f}")
    
    # Format the final answer string as requested
    print("\n" + "#" * 30)
    print("Final Answer for Submission:")
    print(f"Y[{error:.3f}]")
    print("#" * 30)


if __name__ == '__main__':
    titan_calculation()
