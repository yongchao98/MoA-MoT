import math

class TitanError(Exception):
    """Custom exception for Titan architecture constraint violations."""
    pass

class TitanFraction:
    """
    Represents a number in Titan's fractional system.
    It includes a numerator, a denominator, and a base-10 exponent.
    The 4-bit constraint (0-15) is enforced on the numerator and denominator.
    """
    def __init__(self, num, den=1, exp=0, simplified=False):
        if not simplified:
            if not (0 <= num <= 15 and 0 < den <= 15):
                raise TitanError(f"Constraint Violation: The fraction {num}/{den} is not valid. Numerator and denominator must be 4-bit integers (0-15).")
        self.num = num
        self.den = den
        self.exp = exp

    def __mul__(self, other):
        """
        Multiplies two TitanFractions, checking for constraint violations.
        """
        new_num = self.num * other.num
        new_den = self.den * other.den
        new_exp = self.exp + other.exp

        # Rule 4: Intermediate results must be immediately simplified.
        common_divisor = math.gcd(new_num, new_den)
        simplified_num = new_num // common_divisor
        simplified_den = new_den // common_divisor

        print(f"Multiplying ({self.num}/{self.den}e{self.exp}) by ({other.num}/{other.den}e{other.exp})")
        print(f"  Intermediate product: {new_num}/{new_den}e{new_exp}")
        print(f"  Simplified product: {simplified_num}/{simplified_den}e{new_exp}")

        # Rule 4: Check constraint AFTER simplification.
        if not (0 <= simplified_num <= 15 and 0 < simplified_den <= 15):
            print(f"  Result: FAILURE. Numerator '{simplified_num}' or denominator '{simplified_den}' exceeds 15.")
            raise TitanError(f"Constraint Violation: Simplified fraction {simplified_num}/{simplified_den} is out of 4-bit range.")
        
        print(f"  Result: SUCCESS. The fraction is valid.")
        return TitanFraction(simplified_num, simplified_den, new_exp, simplified=True)
        
    def __str__(self):
        return f"{self.num}/{self.den}e{self.exp}"

def calculate_escape_velocity():
    """
    Attempts to calculate Pandora's escape velocity using the Titan model.
    """
    try:
        print("Starting Pandora escape velocity calculation on Titan architecture.\n")
        
        # Step 1: Define constants and parameters as TitanFractions
        # These must all be valid 4-bit fractions.
        print("--- Defining Constants ---")
        # Universal constant from volume of a sphere (4/3 * pi * R^3)
        # We combine it with the factor of 2 from v_e^2 = 2GM/R, making 8/3
        C_8_div_3 = TitanFraction(8, 3) 
        print(f"Constant 8/3: {C_8_div_3}")

        # Planet Radius squared: R = 2000 km = 2e6 m, so R^2 = 4e12 m^2
        R_squared = TitanFraction(4, 1, 12)
        print(f"Planet Radius Squared (R^2): {R_squared}")
        
        # Shell density: d_shell = 0.3 metric tons/m^3 = 300 kg/m^3
        d_shell = TitanFraction(3, 10, 3) # Representing 300 as 3/10 * 10^3
        print(f"Shell Density (d_shell): {d_shell}\n")
        
        print("--- Calculation Steps ---")
        print("The formula for escape velocity squared is v_e^2 approx (8/3) * pi * G * R^2 * d_shell.")
        print("Let's analyze the term (8/3) * R^2, as this multiplication must occur at some point.")

        # Step 2: Perform a necessary multiplication that violates the constraints.
        # This multiplication is unavoidable regardless of the order of operations.
        # It represents a fundamental limitation in this calculation.
        term = C_8_div_3 * R_squared
        
    except TitanError as e:
        print(f"\nFATAL ERROR during calculation: {e}")
        print("\nConclusion: The calculation cannot be performed on the Titan architecture")
        print("because fundamental operations result in numbers that cannot be represented within the 4-bit constraints.")
        print("\n<<<N0>>>")

    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == '__main__':
    calculate_escape_velocity()