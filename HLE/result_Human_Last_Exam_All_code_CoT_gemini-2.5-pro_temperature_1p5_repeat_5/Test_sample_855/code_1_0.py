import math

class TitanComputerError(Exception):
    """Custom exception for Titan computer constraint violations."""
    pass

class TitanNumber:
    """
    Represents a number in the Titan computer's scientific notation format.
    The number is (numerator / denominator) * 10^exponent.
    The numerator and denominator must be 5-bit integers (0-31).
    """
    MAX_VAL = 31

    def __init__(self, num, den=1, exp=0):
        if not (0 <= num <= self.MAX_VAL and 0 < den <= self.MAX_VAL):
            raise TitanComputerError(f"Initialization failed: Numerator ({num}) or denominator ({den}) out of 5-bit range [0, {self.MAX_VAL}]")
        self.num = num
        self.den = den
        self.exp = exp

    def __str__(self):
        return f"({self.num}/{self.den}) * 10^{self.exp}"

    def __mul__(self, other):
        """
        Multiplies two TitanNumbers, respecting the 5-bit constraint for intermediate products.
        Based on the example, we can cancel common factors before multiplying.
        """
        print(f"Attempting to multiply {self} and {other}")

        # Simplify fractions by cancelling common factors before multiplication
        common_factor1 = math.gcd(self.num, other.den)
        common_factor2 = math.gcd(other.num, self.den)

        num1 = self.num // common_factor1
        den2 = other.den // common_factor1
        num2 = other.num // common_factor2
        den1 = self.den // common_factor2
        
        print(f"After cancelling common factors, the operation is ({num1}/{den1}) * ({num2}/{den2})")

        # Calculate the new numerator and denominator from the simplified parts
        new_num = num1 * num2
        new_den = den1 * den2
        
        # Critical Check: The problem states that the result of an operation cannot
        # have a numerator or denominator that exceeds the 5-bit limit.
        print(f"The final equation for the numerator is: {num1} * {num2} = {new_num}")
        print(f"The final equation for the denominator is: {den1} * {den2} = {new_den}")
        
        if new_num > self.MAX_VAL:
            raise TitanComputerError(f"OVERFLOW! Numerator product {new_num} exceeds the 5-bit limit of {self.MAX_VAL}.")
        if new_den > self.MAX_VAL:
             raise TitanComputerError(f"OVERFLOW! Denominator product {new_den} exceeds the 5-bit limit of {self.MAX_VAL}.")

        new_exp = self.exp + other.exp
        return TitanNumber(new_num, new_den, new_exp)

# --- Main Program ---
print("Analyzing feasibility of calculating Landing Force on the Titan Computer.")
print("The force equation is F_rocket = m * (a + g).")
print("We will test the critical sub-calculation of m * a.\n")

try:
    # Represent the probe mass 'm' (50 kg) in Titan format.
    # 50 = 5 * 10^1. The fraction part is 5/1.
    m = TitanNumber(5, 1, 1)
    print(f"Probe mass m = 50 kg is represented as: {m}")

    # The required constant acceleration 'a' is 9 m/s^2.
    # 9 = 9 * 10^0. The fraction part is 9/1.
    a = TitanNumber(9, 1, 0)
    print(f"Required acceleration a = 9 m/s^2 is represented as: {a}\n")

    # Attempt to multiply m and a. This is the core test.
    force_component = m * a
    print(f"\nSuccessfully calculated m*a = {force_component}")
    print("This result suggests the full calculation might be possible.")

except TitanComputerError as e:
    print(f"\n--- CALCULATION FAILED ---")
    print(f"Error Message: {e}")
    print("\nExplanation:")
    print("The multiplication of the fractional parts (5/1 * 9/1) results in a numerator of 45.")
    print("Because 45 is greater than 31, it cannot be represented by a 5-bit integer.")
    print("This violates a fundamental constraint of the Titan architecture.")
    print("\nConclusion: Since a required part of the overall force calculation fails, the entire problem cannot be solved on Titan.")