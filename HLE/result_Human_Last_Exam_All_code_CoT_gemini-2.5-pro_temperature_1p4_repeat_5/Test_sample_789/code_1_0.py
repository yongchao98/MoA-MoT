import math

def main():
    """
    This script calculates the mass of a rock on the Titan 5-bit architecture,
    deriving the result with the smallest possible absolute error.
    """
    MAX_INT = 31

    def gcd(a, b):
        """Helper function to find the greatest common divisor."""
        while b:
            a, b = b, a % b
        return a

    def multiply_fractions(f1_num, f1_den, f2_num, f2_den):
        """
        Multiplies two fractions (f1_num/f1_den) * (f2_num/f2_den)
        respecting the 5-bit constraint by simplifying before multiplying.
        """
        # Simplify across terms to prevent intermediate overflow
        # (f1_num / f2_den) * (f2_num / f1_den)
        common1 = gcd(f1_num, f2_den)
        n1 = f1_num // common1
        d1 = f2_den // common1

        common2 = gcd(f2_num, f1_den)
        n2 = f2_num // common2
        d2 = f1_den // common2

        # Perform the multiplication on simplified terms
        res_num = n1 * n2
        res_den = d1 * d2

        if res_num > MAX_INT or res_den > MAX_INT:
            raise ValueError(f"Overflow: ({f1_num}/{f1_den})*({f2_num}/{f2_den}) -> {res_num}/{res_den}")
        
        return res_num, res_den

    print("Derivation of the rock's mass following Titan 5-bit rules:\n")
    print("Formula: Mass = Density * (4/3) * pi * r^3\n")

    # 1. Represent initial values as fractions
    r_num, r_den = 1, 2  # r = 0.5 cm
    four_thirds_num, four_thirds_den = 4, 3
    pi_num, pi_den = 22, 7  # Most accurate 5-bit approximation for pi
    
    # We strategically approximate density rho=0.9 as 10/11 (~0.909)
    # This enables calculation with the accurate pi value without overflow.
    rho_num, rho_den = 10, 11
    
    print("Step 1: Represent all constants as 5-bit fractions.")
    print(f"  - Radius r = 0.5            -> {r_num}/{r_den}")
    print(f"  - Constant 4/3              -> {four_thirds_num}/{four_thirds_den}")
    print(f"  - Pi (π)                    -> {pi_num}/{pi_den} (best 5-bit approximation)")
    print(f"  - Density (ρ) = 0.9 kg/cm^3 -> {rho_num}/{rho_den} (strategic approximation)\n")

    # 2. Calculate r^3
    r2_num, r2_den = multiply_fractions(r_num, r_den, r_num, r_den)       # (1/2)*(1/2) = 1/4
    r3_num, r3_den = multiply_fractions(r2_num, r2_den, r_num, r_den)   # (1/4)*(1/2) = 1/8
    print(f"Step 2: Calculate r^3 = ({r_num}/{r_den})^3 = {r3_num}/{r3_den}\n")
    
    # 3. Calculate Volume = (4/3) * pi * r^3
    # Grouping (4/3) * r^3 first
    v_const_num, v_const_den = multiply_fractions(four_thirds_num, four_thirds_den, r3_num, r3_den) # (4/3)*(1/8)=1/6
    
    # Then multiply by pi
    v_num, v_den = multiply_fractions(v_const_num, v_const_den, pi_num, pi_den) # (1/6)*(22/7)=11/21
    print(f"Step 3: Calculate Volume = ( (4/3) * (1/8) ) * (22/7) = (1/6) * (22/7) = {v_num}/{v_den}\n")

    # 4. Calculate Mass = Density * Volume
    mass_num, mass_den = multiply_fractions(rho_num, rho_den, v_num, v_den) # (10/11)*(11/21)=10/21
    print(f"Step 4: Calculate Mass = ({rho_num}/{rho_den}) * ({v_num}/{v_den}) = {mass_num}/{mass_den}\n")

    # 5. Final Equation and Error Calculation
    print("--- Summary ---")
    print("The full calculation is performed as follows:")
    final_eq_str = (
        f"Mass = (Density) * (4/3) * (pi) * (r^3)\n"
        f"     = ({rho_num}/{rho_den}) * ({four_thirds_num}/{four_thirds_den}) * ({pi_num}/{pi_den}) * ({r3_num}/{r3_den})\n"
        f"     = {mass_num}/{mass_den}"
    )
    print(final_eq_str)
    
    # Calculate true mass and error
    true_mass = 0.9 * (4.0/3.0) * math.pi * (0.5**3)
    calculated_mass = mass_num / mass_den
    absolute_error = abs(calculated_mass - true_mass)
    
    print(f"\nCalculated Mass = {calculated_mass:.5f} kg")
    print(f"True Mass       = {true_mass:.5f} kg")
    print(f"Absolute Error  = {absolute_error:.5f}")

    # Answer is the error rounded to 0.001
    final_answer = round(absolute_error, 3)
    
    print(f"\nThe smallest absolute error, rounded to 0.001, is: {final_answer}")
    
    print(f"<<<{final_answer}>>>")

if __name__ == "__main__":
    main()