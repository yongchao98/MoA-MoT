import math

def solve_titan_mass_problem():
    """
    This function simulates the calculation of the rock's mass using
    Titan's 5-bit fractional arithmetic rules to find the smallest possible error.
    """
    
    # --- Define a helper function for Titan's multiplication rules ---
    def titan_multiply(f1_num, f1_den, f2_num, f2_den):
        """
        Multiplies two fractions, simplifying before multiplication to keep
        intermediate values small, as required by the Titan architecture.
        Returns the resulting (numerator, denominator) or None if it overflows.
        """
        # Simplify f1_num with f2_den
        common1 = math.gcd(f1_num, f2_den)
        num1 = f1_num // common1
        den2 = f2_den // common1

        # Simplify f2_num with f1_den
        common2 = math.gcd(f2_num, f1_den)
        num2 = f2_num // common2
        den1 = f1_den // common2

        # Multiply the simplified parts
        res_num = num1 * num2
        res_den = den1 * den2

        # On Titan, an overflow would occur if results exceed 5-bit limits
        if res_num > 31 or res_den > 31:
            return None
        
        return res_num, res_den

    # --- Step 1: State initial values as 5-bit fractions ---
    rho_num, rho_den = 9, 10  # density = 0.9
    r_num, r_den = 1, 2      # radius = 0.5
    four_thirds_num, four_thirds_den = 4, 3

    print("Derivation of the mass calculation:")
    print(f"mass = (density) * (4/3) * pi * (radius)^3")
    print(f"mass = ({rho_num}/{rho_den}) * ({four_thirds_num}/{four_thirds_den}) * pi * ({r_num}/{r_den})^3")
    
    # --- Step 2: Calculate r^3 ---
    # (1/2)^3 = 1/8. This is a valid 5-bit fraction.
    r3_num, r3_den = 1, 8
    print(f"First, r^3 = ({r_num}/{r_den})^3 = {r3_num}/{r3_den}")

    # --- Step 3: Combine constant terms to find the coefficient of pi ---
    # Term 1: (9/10) * (4/3) = 36/30. This overflows.
    # We must simplify during multiplication.
    # (9/3) * (4/10) = 3 * (2/5) = 6/5
    temp_num, temp_den = titan_multiply(rho_num, rho_den, four_thirds_num, four_thirds_den)
    print(f"Combine constants: ({rho_num}/{rho_den}) * ({four_thirds_num}/{four_thirds_den}) = {temp_num}/{temp_den}")

    # Term 2: (6/5) * (1/8) = 6/40. This overflows.
    # Simplify during multiplication: (6/2) / (40/2) -> 3/20
    pi_coeff_num, pi_coeff_den = titan_multiply(temp_num, temp_den, r3_num, r3_den)
    print(f"Then multiply by r^3: ({temp_num}/{temp_den}) * ({r3_num}/{r3_den}) = {pi_coeff_num}/{pi_coeff_den}")
    print(f"The expression is now: mass = ({pi_coeff_num}/{pi_coeff_den}) * pi")

    # --- Step 4: Choose pi approximation and calculate final mass ---
    # We need to multiply (3/20) by pi. To avoid overflow, we need an
    # approximation of pi like c/d where c is divisible by a factor of 20
    # and/or d is divisible by 3.
    # pi approx 28/9 is a good candidate (28 is divisible by 4, 9 by 3).
    pi_approx_num, pi_approx_den = 28, 9
    print(f"Approximating pi as {pi_approx_num}/{pi_approx_den} to enable simplification:")
    
    final_num, final_den = titan_multiply(pi_coeff_num, pi_coeff_den, pi_approx_num, pi_approx_den)
    print(f"Final multiplication: ({pi_coeff_num}/{pi_coeff_den}) * ({pi_approx_num}/{pi_approx_den}) simplifies to {final_num}/{final_den}")
    print("-" * 30)

    # --- Step 5: Calculate the error ---
    true_mass = 0.9 * (4.0/3.0) * math.pi * (0.5**3)
    calculated_mass = final_num / final_den
    absolute_error = abs(true_mass - calculated_mass)
    
    print("Final Result:")
    print(f"The final equation for the mass is: {final_num} / {final_den}")
    print(f"The smallest absolute error found is {absolute_error:.5f}")
    
    # The final answer is the error rounded to 0.001
    final_answer = round(absolute_error, 3)
    print(f"Rounded to 0.001, the error is: {final_answer}")
    
    return final_answer

# Execute the calculation and print the final answer in the required format
result = solve_titan_mass_problem()
print(f"<<<{result}>>>")
