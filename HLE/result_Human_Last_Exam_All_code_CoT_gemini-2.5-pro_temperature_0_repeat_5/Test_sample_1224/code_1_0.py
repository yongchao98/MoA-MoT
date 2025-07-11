import math

# This class simulates fractional arithmetic on the Titan 4-bit architecture.
class TitanFraction:
    """
    A class to represent a number in Titan's fractional system.
    It stores a numerator, a denominator, and a base-10 exponent.
    It enforces the 4-bit constraint (0-15) for numerators and denominators.
    """
    def __init__(self, num, den, exp=0):
        if not (0 <= num <= 15 and 0 <= den <= 15):
            raise ValueError(f"Numerator ({num}) or denominator ({den}) exceeds 4-bit limit (0-15).")
        if den == 0:
            raise ValueError("Denominator cannot be zero.")
        
        self.num = num
        self.den = den
        self.exp = exp
        self.simplify()

    def simplify(self):
        """Reduces the fraction by dividing num and den by their GCD."""
        if self.num == 0:
            self.den = 1
            return
        common_divisor = math.gcd(self.num, self.den)
        self.num //= common_divisor
        self.den //= common_divisor

    def __repr__(self):
        """String representation of the fraction."""
        if self.exp == 0:
            return f"{self.num}/{self.den}"
        return f"{self.num}/{self.den} * 10^{self.exp}"
        
    def to_decimal(self):
        return (self.num / self.den) * (10**self.exp)

    def __mul__(self, other):
        """Multiplies two TitanFractions."""
        num = self.num * other.num
        den = self.den * other.den
        exp = self.exp + other.exp
        # This operation would fail if num or den > 15, so we must simplify before creating the new object.
        # For this simulation, we create it and let the __init__ check fail.
        return TitanFraction(num, den, exp)

    def __div__(self, other):
        """Divides two TitanFractions."""
        num = self.num * other.den
        den = self.den * other.num
        exp = self.exp - other.exp
        return TitanFraction(num, den, exp)

    def __add__(self, other):
        """Adds two TitanFractions."""
        # To add, exponents must be the same.
        if self.exp != other.exp:
            # This simulation does not implement exponent alignment for brevity,
            # as the calculation will fail on a more fundamental constraint.
            raise ValueError("Cannot add fractions with different exponents in this simplified simulation.")
        
        num = self.num * other.den + other.num * self.den
        den = self.den * other.den
        exp = self.exp
        return TitanFraction(num, den, exp)

def run_titan_calculation():
    """
    Simulates the calculation of the landing time on Pandora using the Titan architecture.
    """
    try:
        print("--- Titan Feasibility Calculation ---")
        print("Step 1: Define constants using 4-bit fractional approximations.")
        
        # g = 4/3 * pi * G * r * rho. We simplify pi=3, so g = 4 * (G*rho) * r
        # G*rho is approximated as 2e-8
        G_rho = TitanFraction(2, 1, -8)
        # r is approximated as 2000 km = 2e6 m
        r = TitanFraction(2, 1, 6)
        # h = 5000 m
        h = TitanFraction(5, 1, 3)
        
        four = TitanFraction(4, 1)
        two = TitanFraction(2, 1)
        
        print(f"Approximated G*rho = {G_rho} s^-2")
        print(f"Approximated radius r = {r} m")
        print(f"Fall height h = {h} m")
        print("\nStep 2: Calculate gravitational acceleration 'g'.")
        
        # g = 4 * G_rho * r
        g_intermediate = four * G_rho
        print(f"  4 * (G*rho) = {four} * {G_rho} = {g_intermediate}")
        g = g_intermediate * r
        print(f"  g = {g_intermediate} * {r} = {g}")
        print(f"Result: g = {g.to_decimal()} m/s^2. Calculation successful.\n")

        print("Step 3: Calculate t^2 = 2*h/g.")
        two_h = two * h
        print(f"  2 * h = {two} * {h} = {two_h}")
        t_squared = two_h / g
        print(f"  t^2 = {two_h} / {g} = {t_squared}")
        print(f"Result: t^2 = {t_squared.to_decimal()}. Calculation successful.\n")

        print("Step 4: Attempt to calculate t = sqrt(t^2) using Newton-Raphson method.")
        print("The method is x_next = 0.5 * (x_prev + t^2 / x_prev).")
        
        # Initial guess for sqrt(62500) is 200.
        x0 = TitanFraction(2, 1, 2) # 200
        print(f"Initial guess x0 = {x0}")
        
        # Calculate t^2 / x0
        t_sq_div_x0 = t_squared / x0
        print(f"  t^2 / x0 = {t_squared} / {x0} = {t_sq_div_x0}")
        
        # Now, add x0 + (t^2 / x0). This requires common exponents.
        # x0 = 2e2. t_sq_div_x0 = 5/16 e3 = 50/16 e2 = 25/8 e2.
        # The numerator 25 is > 15, so this representation is invalid.
        print("\n--- FAILURE ---")
        print("To add x0 and (t^2/x0), their exponents must match.")
        print(f"Let's convert {t_sq_div_x0} to have an exponent of 2.")
        print(f"  {t_sq_div_x0} = 5/16 * 10^3 = 50/16 * 10^2")
        print("  Simplifying 50/16 gives 25/8.")
        
        # This is the point of failure.
        invalid_fraction = TitanFraction(25, 8, 2)

    except ValueError as e:
        print(f"An operation resulted in a value that cannot be represented on Titan: {e}")
        print("The calculation cannot be completed due to the 4-bit integer limitation.")
        print("\nConclusion: It is not possible to use Titan to calculate this landing time.")
        print("<<<N0>>>")

run_titan_calculation()