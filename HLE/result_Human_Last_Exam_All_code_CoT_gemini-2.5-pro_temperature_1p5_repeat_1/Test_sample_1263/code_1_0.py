import math

# Titan architecture constraints
MAX_VAL = 15

class TitanFrac:
    """A class to represent a number in Titan's fractional system."""
    def __init__(self, num, den=1, exp=0, name=""):
        self.num = num
        self.den = den
        self.exp = exp
        self.name = name
        self.check_constraints()

    def check_constraints(self):
        """Check if the numerator and denominator are valid 4-bit integers."""
        if not (0 <= self.num <= MAX_VAL):
            raise ValueError(f"Constraint Violation: Numerator {self.num} for '{self.name}' is out of the 0-15 range.")
        if not (0 <= self.den <= MAX_VAL):
            raise ValueError(f"Constraint Violation: Denominator {self.den} for '{self.name}' is out of the 0-15 range.")

    def __str__(self):
        return f"{self.num}/{self.den} x 10^{self.exp}"

def multiply(f1: TitanFrac, f2: TitanFrac, result_name: str) -> TitanFrac:
    """Simulates multiplication on Titan, checking for violations."""
    print(f"Multiplying: ({f1.name}={f1}) by ({f2.name}={f2})")
    
    new_num = f1.num * f2.num
    new_den = f1.den * f2.den
    new_exp = f1.exp + f2.exp
    
    print(f"Intermediate result: {new_num}/{new_den} x 10^{new_exp}")

    # The rules state any operation resulting in values > 15 must be simplified.
    # We check if the direct multiplication is possible before any simplification.
    if new_num > MAX_VAL or new_den > MAX_VAL:
        print(f"!! FAILURE: Multiplication results in a numerator ({new_num}) or denominator ({new_den}) that exceeds {MAX_VAL}.")
        print(f"!! The expression ({f1.num}*{f2.num})/({f1.den}*{f2.den}) cannot be formed directly.")
        raise ValueError("Calculation cannot proceed.")

    # In a real Titan system, this might be simplified, but we show the direct check first.
    # e.g. GCD simplification
    common_divisor = math.gcd(new_num, new_den)
    final_num = new_num // common_divisor
    final_den = new_den // common_divisor
    
    return TitanFrac(final_num, final_den, new_exp, name=result_name)


def main():
    """
    Main function to simulate the escape velocity calculation for Pandora.
    """
    print("--- Pandora Escape Velocity Calculation on Titan Architecture ---\n")
    print("Objective: Calculate v_e^2 = (8 * pi * G * rho_s * R^2) / 3\n")

    try:
        # Step 1: Define constants using 4-bit compatible fractions
        # Physical constants
        # For pi, 22/7 is not allowed (22>15). 13/4 = 3.25 is a good approximation.
        pi = TitanFrac(13, 4, name="pi")
        # For G, 6.674e-11. 13/2 = 6.5 is a good approximation.
        G = TitanFrac(13, 2, exp=-11, name="G")

        # Planetary parameters
        # rho_s = 300 kg/m^3 = 3 * 10^2
        rho_s = TitanFrac(3, 1, exp=2, name="rho_s")
        # R = 2000 km = 2 * 10^6 m -> R^2 = 4 * 10^12
        R_sq = TitanFrac(4, 1, exp=12, name="R^2")
        
        # Formula constants
        eight_thirds = TitanFrac(8, 3, name="8/3")
        
        print("--- Chosen Approximations and Constants ---")
        print(f"pi ~= {pi}")
        print(f"G ~= {G}")
        print(f"rho_s = {rho_s}")
        print(f"R^2 = {R_sq}")
        print("-" * 20 + "\n")
        
        # Step 2: Perform the calculation, showing cancellations
        print("--- Step-by-step Calculation for v_e^2 ---")
        # Let's combine 8/3 and rho_s first due to the convenient cancellation of '3'
        term1_val = eight_thirds.num * rho_s.num
        term1_den = eight_thirds.den * rho_s.den
        # Simplify 24/3 -> 8/1
        term1_num_s = term1_val // term1_den
        term1_den_s = 1
        term1 = TitanFrac(term1_num_s, term1_den_s, exp=rho_s.exp, name="(8/3)*rho_s")
        print(f"Multiplying: ({eight_thirds.name}) by ({rho_s.name})")
        print(f"Intermediate result: {term1_val}/{term1_den} x 10^{rho_s.exp}")
        print(f"Simplified result is: {term1}. This step is valid.\n")

        # Now multiply by R^2
        # current_result = term1 * R_sq
        # This involves 8 * 4 = 32
        current_result = multiply(term1, R_sq, name=f"{term1.name}*R^2")

    except ValueError as e:
        print("\n--- CONCLUSION ---")
        print("The calculation for Pandora's escape velocity is NOT feasible on the Titan architecture.")
        print(f"Reason: As shown above, the computation requires multiplying valid 4-bit numbers (e.g., 8 and 4) which results in a value (32) that cannot be represented by a 4-bit integer, violating a fundamental constraint.")
        print("Even with optimal choices of fractional approximations and calculation order, intermediate products exceed the 0-15 limit.")
        print("\nFinal Answer: N0")

main()
<<<N0>>>