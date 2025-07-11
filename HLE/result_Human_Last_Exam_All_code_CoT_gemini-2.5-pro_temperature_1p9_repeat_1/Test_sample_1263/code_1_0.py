import math

class TitanComputer:
    """
    A simulator for the Titan 4-bit computer architecture.
    It defines rules for fractional arithmetic under Titan's constraints.
    """

    class Fraction:
        def __init__(self, num, den=1, exponent=0):
            if not isinstance(num, int) or not isinstance(den, int):
                raise TypeError("Numerator and denominator must be integers.")
            if den == 0:
                raise ValueError("Denominator cannot be zero.")

            # Auto-simplification as part of representation
            common = math.gcd(num, den)
            self.num = num // common
            self.den = den // common
            self.exponent = exponent
            self.validate()

        def validate(self):
            """Check if the fraction adheres to the 4-bit constraint (0-15)."""
            if not (0 <= abs(self.num) <= 15 and 0 <= abs(self.den) <= 15):
                raise ValueError(f"Constraint Violated: Fraction {self.num}/{self.den} is outside the 4-bit limit.")
            return True

        def __mul__(self, other):
            new_num = self.num * other.num
            new_den = self.den * other.den
            new_exp = self.exponent + other.exponent
            # This represents the operation and immediate simplification
            return TitanComputer.Fraction(new_num, new_den, new_exp)

        def __truediv__(self, other):
            new_num = self.num * other.den
            new_den = self.den * other.num
            new_exp = self.exponent - other.exponent
            return TitanComputer.Fraction(new_num, new_den, new_exp)

        def __str__(self):
            if self.exponent == 0:
                return f"{self.num}/{self.den}"
            else:
                return f"{self.num}/{self.den} * 10^{self.exponent}"

def calculate_escape_velocity():
    """
    Attempts to calculate Pandora's escape velocity using the Titan simulator.
    """
    print("Task: Calculate Pandora's escape velocity v_e = sqrt(2 * G * M / R)")
    print("-" * 20)
    print("Step 1: Define constants as Titan Fractions.")

    try:
        # Define constants using best 4-bit approximations
        # Note: Exponent for M is chosen to make the final exponent of v_e^2 even.
        const_2 = TitanComputer.Fraction(2, 1)
        G = TitanComputer.Fraction(13, 2, exponent=-11)  # G ≈ 6.5e-11
        M = TitanComputer.Fraction(10, 1, exponent=21)  # M ≈ 10e21
        R = TitanComputer.Fraction(2, 1, exponent=6)   # R = 2e6

        print(f"Constant 2 = {const_2}")
        print(f"Approximation for G = {G}")
        print(f"Approximation for M = {M}")
        print(f"Value for R = {R}")
        print("-" * 20)
        
        # v_e^2 = 2 * G * M / R
        print("Step 2: Calculate v_e^2 term by term, mimicking Titan registers.")
        print(f"Equation for v_e^2 = ({const_2}) * ({G}) * ({M}) / ({R})")

        # AX = 2 * G
        print("\nExecuting: MUL AX, G where AX = 2/1")
        print(f"AX = {const_2} * {G}")
        ax = const_2 * G
        print(f"Intermediate result for 2*G: {ax}. This is valid.")

        # AX = (2*G) * M
        print("\nExecuting: MUL AX, M where AX = 13/1 * 10^-11")
        print(f"AX = {ax} * {M}")
        ax = ax * M
        print(f"Intermediate result for 2*G*M: {ax}. This is valid after simplification (130/1 -> 130/1).")
        # The above print will not be reached if validation is strict.
        # Let's show the failure point explicitly.
        
    except ValueError as e:
        print("\n--- CALCULATION FAILED ---")
        print("Failure occurred during the multiplication of (2*G) and M.")
        
        # Let's manually show the operation that fails
        # term1 = 2 * G => (2/1 * 13/2) = 26/2 which simplifies to 13/1 (Valid)
        term1_num = 2 * 13
        term1_den = 1 * 2
        
        # term2 = (2*G) * M => (13/1 * 10/1) = 130/1
        term2_num = (term1_num // math.gcd(term1_num, term1_den)) * 10
        term2_den = (term1_den // math.gcd(term1_num, term1_den)) * 1

        print(f"The operation is (13/1) * (10/1).")
        print(f"The resulting fraction is {term2_num}/{term2_den}.")
        print(f"Reason for failure: The numerator '{term2_num}' is greater than 15.")
        print("The 4-bit constraint is violated, and the number cannot be represented on Titan.")
        print("\nSince a critical part of the calculation is impossible, Titan cannot be used.")
        print("Final Answer: N0")
        
        # Final required output format
        print("\n<<<N0>>>")
        return

if __name__ == '__main__':
    calculate_escape_velocity()