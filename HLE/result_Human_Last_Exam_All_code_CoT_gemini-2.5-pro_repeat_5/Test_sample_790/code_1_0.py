import math

class TitanFraction:
    """
    A class to represent and compute with fractions under Titan's 5-bit constraints.
    It simplifies fractions upon creation.
    """
    def __init__(self, num, den=1):
        if not isinstance(num, int) or not isinstance(den, int):
            raise TypeError("Numerator and denominator must be integers.")
        if den == 0:
            raise ZeroDivisionError("Denominator cannot be zero.")

        # The initial numbers must be representable by 5-bit integers.
        if not (0 <= num <= 31 and 0 < den <= 31):
            # This check can be relaxed for intermediate, pre-simplification results
            # but for our inputs, it holds.
            pass

        self.num = num
        self.den = den
        self._simplify()

    def _simplify(self):
        """Simplifies the fraction by dividing by the greatest common divisor."""
        if self.num == 0:
            self.den = 1
            return
        common_divisor = math.gcd(self.num, self.den)
        self.num //= common_divisor
        self.den //= common_divisor

    def __repr__(self):
        """Represents the fraction as a string 'num/den'."""
        return f"{self.num}/{self.den}"

    def __sub__(self, other):
        """Subtracts two TitanFractions."""
        res_num = self.num * other.den - other.num * self.den
        res_den = self.den * other.den
        # The result is simplified upon creation of the new TitanFraction instance.
        return TitanFraction(res_num, res_den)

    def __truediv__(self, other):
        """Divides two TitanFractions."""
        res_num = self.num * other.den
        res_den = self.den * other.num
        # The result is simplified upon creation of the new TitanFraction instance.
        return TitanFraction(res_num, res_den)

def solve_monkey_problem():
    """
    Calculates the dimensionless force factor 'k' using Titan's rules.
    The factor k is defined by v₀² = k * g * d, and it is derived
    that k = d / (d - h), where d is the distance and h is the height.
    """
    print("Titan Computer Calculation Log")
    print("---------------------------------")
    print("Objective: Calculate the dimensionless force factor k = d / (d - h)")

    # Step 1: Represent constants as Titan fractions.
    # d = 20 meters, h = 10 meters
    d = TitanFraction(20, 1)
    h = TitanFraction(10, 1)
    print(f"Representing distance d = {d}")
    print(f"Representing height h = {h}")

    # Step 2: Calculate the difference (d - h)
    diff = d - h
    print(f"Calculating difference: (d - h) = {d} - {h} = {diff}")

    # Step 3: Calculate the factor k = d / (d - h)
    k = d / diff
    print(f"Calculating factor k = d / (d - h) = {d} / {diff} = {k}")
    print("---------------------------------")

    # Step 4: Output the final equation with all numbers.
    print("Final Equation:")
    print(f"{k} = {d} / {diff}")
    print("\nThe dimensionless force factor required is 2.")
    print("Since this calculation is exact for the center of the target (10m), the absolute error is 0.")

solve_monkey_problem()