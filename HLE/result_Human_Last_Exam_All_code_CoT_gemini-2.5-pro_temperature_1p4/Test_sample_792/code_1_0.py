import math

# Titan Computer Simulation for the Curious Monkey Problem

def get_gcd(a, b):
    """Helper function to calculate the greatest common divisor."""
    return math.gcd(a, b)

class TitanFraction:
    """
    Represents a fraction in the Titan 5-bit architecture.
    Numerators and denominators must be between 0 and 31.
    """
    MAX_INT = 31

    def __init__(self, numerator, denominator):
        if not (0 <= numerator <= self.MAX_INT and 0 <= denominator <= self.MAX_INT):
            raise ValueError(
                f"Numerator ({numerator}) or denominator ({denominator}) "
                f"out of 5-bit range (0-{self.MAX_INT})."
            )
        if denominator == 0:
            raise ZeroDivisionError("Denominator cannot be zero.")
        
        # Simplify the fraction upon creation
        common = get_gcd(numerator, denominator)
        self.num = numerator // common
        self.den = denominator // common

    def __mul__(self, other):
        """
        Multiplies two TitanFractions, simplifying before multiplication
        to prevent intermediate overflow.
        a/b * c/d = (a/gcd(a,d) * c/gcd(c,b)) / (b/gcd(c,b) * d/gcd(a,d))
        """
        g1 = get_gcd(self.num, other.den)
        g2 = get_gcd(other.num, self.den)

        num1 = self.num // g1
        den2 = other.den // g1
        
        num2 = other.num // g2
        den1 = self.den // g2

        new_num = num1 * num2
        new_den = den1 * den2
        
        return TitanFraction(new_num, new_den)

    def __truediv__(self, other):
        """Divides two TitanFractions by multiplying by the reciprocal."""
        reciprocal = TitanFraction(other.den, other.num)
        return self * reciprocal

    def __str__(self):
        return f"{self.num}/{self.den}"

    def to_decimal(self):
        return self.num / self.den

def solve_monkey_problem():
    """
    Solves the problem using Titan's computational rules.
    """
    print("--- Titan Superconducting Computer Calculation ---")
    print("Goal: Calculate the force F needed to throw a rock at a lion.")
    
    # Physics formula: F = (x * m * g) / h
    # With m = (3/20)*pi, x = 20, h = 10, this simplifies to:
    # F = (20/10) * (3/20 * pi) * g = (3/10) * pi * g
    print("\nDerived formula for Titan: F = (3/10) * pi * g")

    # Step 1: Define constants as TitanFractions
    # We choose approximations that allow computation within 5-bit constraints.
    # pi is approximated as 28/9 (~3.111)
    # g is approximated as 10/1 (10.0 m/s^2)
    # The rock's mass constant is 3/20
    m_base = TitanFraction(3, 20)
    pi_approx = TitanFraction(28, 9)
    g_approx = TitanFraction(10, 1)

    print("\nChosen Fractional Approximations:")
    print(f"  - Rock mass constant: {m_base}")
    print(f"  - Pi approximation (pi): {pi_approx}")
    print(f"  - Gravity approximation (g): {g_approx}")

    # Step 2: Perform the calculation step-by-step
    print("\nCalculation Steps:")
    
    # First, calculate (3/10) * pi
    term1 = m_base # This represents the (3/10) factor
    print(f"Step 1: Start with the mass and pi factor: {term1}")
    
    term2 = term1 * pi_approx
    print(f"Step 2: Multiply by pi ({pi_approx}): {term1} * {pi_approx} = {term2}")

    # Then, multiply the result by g
    final_force = term2 * g_approx
    print(f"Step 3: Multiply by g ({g_approx}): {term2} * {g_approx} = {final_force}")
    
    # Final Result
    print("\n--- Final Result ---")
    print(f"The calculated force F required is the fraction:")
    print(f"{final_force.num} / {final_force.den}")
    print(f"This is approximately {final_force.to_decimal():.4f} N.")

    # Error calculation
    print("\n--- Error Analysis ---")
    # F_true = (x/h) * m * g = (20/10) * (3/20 * math.pi) * 9.8
    f_true = (20.0/10.0) * (3.0/20.0 * math.pi) * 9.8
    f_titan = final_force.to_decimal()
    absolute_error = abs(f_true - f_titan)

    print(f"True value of Force (using standard precision): {f_true:.4f} N")
    print(f"Absolute error |true - titan| = |{f_true:.4f} - {f_titan:.4f}| = {absolute_error:.4f}")
    
    # Final Answer Format
    # The calculation is possible, so the answer is Y[e]
    # where e is the smallest absolute error.
    answer = f"Y[{absolute_error:.4f}]"
    print("\nFinal Answer Payload:")
    print(f"<<<{answer}>>>")


solve_monkey_problem()