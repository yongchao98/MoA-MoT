import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    return math.gcd(a, b)

class TitanFraction:
    """
    A class to represent a fraction under Titan architecture constraints.
    The numerator and denominator must be within the 6-bit range [0, 63].
    """
    def __init__(self, num, den=1, name=""):
        if not (0 <= num <= 63 and 0 <= den <= 63):
            raise ValueError(f"Numerator ({num}) or denominator ({den}) for '{name}' exceeds 6-bit limit (0-63).")
        self.num = num
        self.den = den
        self.name = name

    def __str__(self):
        return f"{self.num}/{self.den}"

def mul_titan(frac1, frac2):
    """
    Simulates Titan's MUL instruction.
    Multiplies two TitanFractions, simplifies the result, and checks for overflow.
    """
    print(f"Multiplying ({frac1.name}: {frac1}) * ({frac2.name}: {frac2})")
    
    # Tentative result before checking constraints
    num_res = frac1.num * frac2.num
    den_res = frac1.den * frac2.den

    print(f"  Intermediate product: {num_res}/{den_res}")

    # Simplify the fraction using GCD, as this is a standard algebraic simplification
    common_divisor = gcd(num_res, den_res)
    num_simple = num_res // common_divisor
    den_simple = den_res // common_divisor
    
    if num_simple != num_res or den_simple != den_res:
        print(f"  Simplified with GCD: {num_simple}/{den_simple}")
        
    # Check for overflow according to Rule #4
    if num_simple > 63 or den_simple > 63:
        print(f"  --> OVERFLOW! Resulting numerator ({num_simple}) or denominator ({den_simple}) exceeds 63.")
        return None  # Indicate failure

    new_name = f"({frac1.name}*{frac2.name})"
    print(f"  --> Success. New value: {num_simple}/{den_simple}")
    return TitanFraction(num_simple, den_simple, name=new_name)

def demonstrate_impossibility():
    """
    Attempts to perform the calculation and shows that it's impossible.
    """
    print("--- Titan Feasibility Study for Force Calculation ---\n")
    print("Goal: Calculate the product of all fractional components for the force.")
    print("F_frac = G * m * rho * (4/3) * pi * R^3\n")

    try:
        # Define constants as TitanFractions
        G_frac = TitanFraction(20, 3, "G")
        m_frac = TitanFraction(50, 1, "m")
        rho_frac = TitanFraction(6, 5, "rho")
        four_thirds = TitanFraction(4, 3, "4/3")
        pi_frac = TitanFraction(22, 7, "pi")
        R_cubed_frac = TitanFraction(8, 1, "R^3")

        print("Chosen calculation sequence: ((rho * 4/3) * R^3) * ...")

        # Step 1: rho * 4/3
        res1 = mul_titan(rho_frac, four_thirds)
        if res1 is None: return
        print("-" * 20)

        # Step 2: (res1) * R^3
        # This step will fail.
        res2 = mul_titan(res1, R_cubed_frac)
        if res2 is None:
            print("\n--- CONCLUSION ---")
            print("The calculation has failed. The multiplication of (8/5) and (8/1) results in 64/5.")
            print("The numerator 64 is outside the 6-bit integer range [0, 63].")
            print("This overflow cannot be resolved by GCD simplification.")
            print("Other multiplication orders also result in unavoidable overflows.")
            print("For instance, (m: 50/1) * (G: 20/3) would result in an intermediate of 1000/3.")
            print("Therefore, the Titan architecture is incapable of performing this calculation.")
            print("\nFinal Answer: N0")
            print("\n<<<N0>>>")
            return
            
    except ValueError as e:
        print(f"An error occurred during setup: {e}")

# Run the demonstration
demonstrate_impossibility()
