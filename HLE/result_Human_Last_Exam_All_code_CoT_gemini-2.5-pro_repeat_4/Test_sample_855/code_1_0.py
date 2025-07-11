import math

# Helper class to simulate Titan's 5-bit fractional representation and arithmetic.
# It checks for overflow (numerator or denominator > 31) during operations.
class TitanFraction:
    def __init__(self, num, den, exp=0):
        if not all(isinstance(i, int) for i in [num, den, exp]):
            raise TypeError("Numerator, denominator, and exponent must be integers.")
        # In a real Titan, these checks would happen at every step.
        # We simplify by checking at creation and during operations.
        if abs(num) > 31 or abs(den) > 31:
            raise ValueError(f"Numerator/Denominator out of 5-bit range: {num}/{den}")
        if den == 0:
            raise ZeroDivisionError("Denominator cannot be zero.")
        
        self.num = num
        self.den = den
        self.exp = exp

    def __mul__(self, other):
        # Multiply mantissas and add exponents
        new_num = self.num * other.num
        new_den = self.den * other.den
        new_exp = self.exp + other.exp

        if abs(new_num) > 31 or abs(new_den) > 31:
            # This operation would be invalid on Titan without simplification.
            # We will find a path that avoids this error.
            raise ValueError(f"Overflow in multiplication: {new_num}/{new_den}")

        return TitanFraction(new_num, new_den, new_exp)

    def __truediv__(self, other):
        # Invert and multiply
        inv_other = TitanFraction(other.den, other.num, -other.exp)
        return self.__mul__(inv_other)

    def to_float(self):
        return (self.num / self.den) * (10 ** self.exp)

    def __str__(self):
        if self.exp != 0:
            return f"({self.num}/{self.den} * 10^{self.exp})"
        return f"({self.num}/{self.den})"

def solve():
    """
    Solves the Pandora landing problem using simulated Titan computer logic.
    """
    print("--- Titan Landing Computer ---")
    print("\nStep 1: Calculate required deceleration a_net = v0^2 / (2*h)")

    # Define initial parameters as Titan fractions
    # v0 = 300 m/s = 3 * 10^2
    v0 = TitanFraction(3, 1, 2)
    print(f"Initial velocity v0 = {v0.to_float()} m/s, represented as {v0}")
    
    # h = 5000 m = 5 * 10^3
    h = TitanFraction(5, 1, 3)
    print(f"Altitude h = {h.to_float()} m, represented as {h}")

    # Calculate v0^2
    v0_sq = TitanFraction(9, 1, 4) # (3*3, 1*1, 2+2)
    print(f"v0^2 = {v0_sq.to_float()}, represented as {v0_sq}")

    # Calculate 2*h
    two_h = TitanFraction(10, 1, 3) # (2*5, 1*1, 0+3)
    # Simplify to keep mantissa small: 10*10^3 -> 1*10^4
    two_h_simplified = TitanFraction(1, 1, 4)
    print(f"2*h = {two_h_simplified.to_float()}, represented as {two_h_simplified}")

    # Calculate a_net = v0^2 / (2*h)
    a_net_num = v0_sq.num * two_h_simplified.den
    a_net_den = v0_sq.den * two_h_simplified.num
    a_net_exp = v0_sq.exp - two_h_simplified.exp
    a_net = TitanFraction(a_net_num, a_net_den, a_net_exp)
    
    print(f"Calculated deceleration a_net = {a_net.to_float()} m/s^2, represented as {a_net}")
    
    print("\nStep 2: Analyze gravitational acceleration g")
    print("Calculating g involves many steps and large intermediate numbers (e.g., 4/3 * pi * R^3).")
    print("A high-precision calculation shows g is approximately 0.17 m/s^2.")
    print(f"Comparing g (~0.17) to a_net (9.0), we see that g is only about 2% of a_net.")
    print("According to Titan Rule 5 (Result Reduction), we can eliminate this negligible term.")
    print("This simplification introduces an error equal to the dropped gravitational force.")
    
    print("\nStep 3: Calculate the simplified landing force F = m * a_net")
    # m_probe = 50 kg
    # Standard representation 5 * 10^1 leads to overflow: (5/1)*(9/1) -> 45/1.
    # We must find a valid representation. Let's try m = 0.5 * 10^2
    m_probe = TitanFraction(1, 2, 2)
    print(f"Probe mass m = {m_probe.to_float()} kg, represented as {m_probe}")
    
    # F = m_probe * a_net
    F_num = m_probe.num * a_net.num
    F_den = m_probe.den * a_net.den
    F_exp = m_probe.exp + a_net.exp
    F_rocket = TitanFraction(F_num, F_den, F_exp)
    
    print("\nFinal Calculation:")
    print(f"F_rocket = m_probe * a_net")
    # The final equation with each number outputted
    print(f"F_rocket = {str(m_probe)} * {str(a_net)} = {str(F_rocket)}")
    final_force_value = F_rocket.to_float()
    print(f"Resulting Force = {final_force_value} N")

    print("\nStep 4: Calculate the absolute error 'e'")
    # Error is the force of gravity we ignored.
    # Use high-precision values for error calculation.
    G = 6.674e-11
    R_p = 2e6 # equatorial radius
    r_c = 1e5 # core radius
    rho_s = 300 # shell density
    rho_c = 1200 # core density
    
    # M_total = M_core + M_shell
    V_core = (4/3) * math.pi * r_c**3
    V_shell = (4/3) * math.pi * (R_p**3 - r_c**3)
    M_total = (rho_c * V_core) + (rho_s * V_shell)
    
    g_true = (G * M_total) / R_p**2
    e = 50 * g_true # The force we ignored
    
    print(f"The true value of g is approx {g_true:.4f} m/s^2.")
    print(f"The absolute error 'e' is the ignored gravitational force F_g = m * g.")
    print(f"e = 50 kg * {g_true:.4f} m/s^2 = {e:.3f} N")

    # Final answer format
    print(f"\nFinal Answer String: Y[{e:.3f}]")
    print(f"\n<<<Y[{e:.3f}]>>>")

solve()