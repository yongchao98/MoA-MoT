import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    return math.gcd(a, b)

def simplify_fraction(n, d):
    """Simplifies a fraction by dividing numerator and denominator by their GCD."""
    if d == 0:
        raise ValueError("Denominator cannot be zero.")
    common = gcd(n, d)
    return n // common, d // common

def check_5bit(num):
    """Checks if a number fits in a 5-bit unsigned integer (0-31)."""
    return 0 <= num <= 31

class TitanComputer:
    """
    A class to simulate the Titan Computer's arithmetic logic.
    """
    def __init__(self):
        self.error_log = []

    def multiply(self, f1, f2):
        """
        Multiplies two fractions f1=(n1, d1) and f2=(n2, d2)
        respecting Titan's simplification and 5-bit constraints.
        """
        n1, d1 = f1
        n2, d2 = f2

        # --- Simplification Step ---
        # Simplify n1 and d2
        common1 = gcd(n1, d2)
        n1_s, d2_s = n1 // common1, d2 // common1
        # Simplify n2 and d1
        common2 = gcd(n2, d1)
        n2_s, d1_s = n2 // common2, d1 // common2
        
        # --- Constraint Check Step ---
        new_n = n1_s * n2_s
        new_d = d1_s * d2_s
        
        print(f"Attempting to multiply ({n1}/{d1}) * ({n2}/{d2})")
        print(f"   After cross-simplification: ({n1_s}/{d1_s}) * ({n2_s}/{d2_s})")
        print(f"   Resulting numerator: {new_n}")
        print(f"   Resulting denominator: {new_d}")


        if not check_5bit(new_n) or not check_5bit(new_d):
            self.error_log.append(
                f"INVALID OPERATION: Result ({new_n}/{new_d}) exceeds 5-bit limit."
            )
            return None
        
        print(f"   SUCCESS: Result is ({new_n}/{new_d})\n")
        return (new_n, new_d)

def solve():
    """
    Main function to perform the calculation based on Titan's architecture.
    """
    titan = TitanComputer()
    
    print("Step 1: Calculate the mass of Pandora's shell (M_shell).")
    print("Formula: M_shell = (4/3) * pi * R^3 * rho_shell")
    print("We will calculate the fractional part of the equation first.\n")

    # Define constants as Titan fractions
    # Using approximations that are valid within 5-bit constraints.
    # pi ≈ 25/8 = 3.125
    # Shell outer radius R = 2000 km = 2e6 m, R^3 = 8e18 m^3
    # Shell density rho_shell = 300 kg/m^3 = 3 * 10^2 kg/m^3
    
    four_thirds = (4, 3)
    pi = (25, 8)
    R_cubed_factor = (8, 1) # Just the fractional part of the number
    rho_shell_factor = (3, 1) # Just the fractional part of the number

    print(f"The equation for the fractional part is: ({four_thirds[0]}/{four_thirds[1]}) * ({pi[0]}/{pi[1]}) * ({R_cubed_factor[0]}/{R_cubed_factor[1]}) * ({rho_shell_factor[0]}/{rho_shell_factor[1]})")
    print("-" * 20)
    
    # Let's group operations strategically to see if a valid path exists.
    # Path 1: ( (4/3) * (3/1) ) * ( (25/8) * (8/1) )
    print("Chosen operational order: [ (4/3) * (3/1) ] * [ (25/8) * (8/1) ]\n")
    
    # First multiplication
    term1 = titan.multiply(four_thirds, rho_shell_factor)
    if not term1:
        print("Calculation failed at first multiplication.")
        print("\nConclusion: Pandora's mass cannot be calculated. Therefore, the landing force cannot be determined.")
        print("<<<N0>>>")
        return

    # Second multiplication
    term2 = titan.multiply(pi, R_cubed_factor)
    if not term2:
        print("Calculation failed at second multiplication.")
        print("\nConclusion: Pandora's mass cannot be calculated. Therefore, the landing force cannot be determined.")
        print("<<<N0>>>")
        return
        
    # Final multiplication
    print("Final step: Multiply the results of the previous operations.")
    mass_factor = titan.multiply(term1, term2)
    
    if not mass_factor:
        print("\n" + titan.error_log[-1])
        print("\nConclusion: Pandora's mass cannot be calculated due to the 5-bit constraint on multiplication results.")
        print("Since the planet's mass is a required component for calculating gravity, the final landing force cannot be determined under Titan's computational rules.")
        print("<<<N0>>>")
    else:
        # This part should not be reached based on the logic.
        print(f"Successfully calculated mass factor: {mass_factor}")
        print("...continuing calculation... (this indicates an unexpected success)")

solve()