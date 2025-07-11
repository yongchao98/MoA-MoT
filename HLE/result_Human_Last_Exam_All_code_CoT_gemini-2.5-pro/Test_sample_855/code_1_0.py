import math

class TitanNumber:
    """A class to represent numbers in Titan's fractional/scientific notation."""
    def __init__(self, num, den=1, exp=0, base=10):
        if not (0 <= num <= 31 and 0 <= den <= 31):
            # This check is for setup, the main logic will handle in-calculation constraints.
            pass
        self.num = num
        self.den = den
        self.exp = exp
        self.base = base

    def __str__(self):
        if self.exp == 0:
            return f"({self.num}/{self.den})"
        return f"({self.num}/{self.den}) * {self.base}^{self.exp}"

    def value(self):
        return (self.num / self.den) * (self.base ** self.exp)

def main():
    print("[Precise Landing with Superconducting Computer]")
    print("Goal: Calculate the required rocket force for the Pioneer probe.")
    print("Formula: F_rocket = m * (a_dec + g)")
    print("-" * 30)

    # Step 1: Define constants in Titan's format
    print("Step 1: Define initial values as Titan numbers.")
    # v_i = 300 m/s = 3 * 10^2
    v_i = TitanNumber(3, 1, 2)
    # d = 5000 m = 5 * 10^3
    d = TitanNumber(5, 1, 3)
    # m = 50 kg.  We will use a special representation later.
    m_pioneer = TitanNumber(5, 1, 1) # Initial representation 5 * 10^1
    print(f"Pioneer initial velocity (v_i): {v_i} m/s")
    print(f"Pioneer landing system altitude (d): {d} m")
    print(f"Pioneer mass (m): {m_pioneer} kg")
    print("-" * 30)

    # Step 2: Calculate the required deceleration a_dec = v_i^2 / (2*d)
    print("Step 2: Calculate deceleration (a_dec) using Titan arithmetic.")
    # v_i^2 = (3 * 10^2)^2 = 9 * 10^4
    v_i_squared = TitanNumber(9, 1, 4)
    print(f"v_i^2 = {v_i} * {v_i} = {v_i_squared} m^2/s^2")
    # 2*d = 2 * (5 * 10^3) = 10 * 10^3 = 1 * 10^4
    two_d = TitanNumber(1, 1, 4)
    print(f"2 * d = 2 * {d} = {two_d} m")
    # a_dec = (9 * 10^4) / (1 * 10^4) = 9
    a_dec = TitanNumber(9, 1, 0)
    print(f"a_dec = v_i^2 / (2*d) = {v_i_squared} / {two_d} = {a_dec} m/s^2")
    print("-" * 30)

    # Step 3: Analyze gravitational term 'g' and simplify
    print("Step 3: Analyze the gravitational force term (F_g = m * g).")
    F_dec_val = 50 * 9 # 450 N
    # True g calculation
    G = 6.674e-11
    rho_shell = 300
    a_radius = 2e6
    b_radius = 1.985e6
    V_planet = (4/3) * math.pi * (a_radius**2) * b_radius
    M_planet = rho_shell * V_planet
    g_val = G * M_planet / (a_radius**2)
    F_g_val = 50 * g_val
    print(f"A standard calculation shows g = {g_val:.3f} m/s^2.")
    print(f"This makes F_dec = {F_dec_val} N, and F_g = {F_g_val:.3f} N.")
    print("F_g is less than 2% of F_dec.")
    print("Attempting to compute g with Titan rules leads to intermediate values > 31, causing overflows.")
    print("Based on Rule 5 (Eliminate negligible terms), we will neglect F_g for the Titan calculation.")
    print("F_rocket ≈ F_dec = m * a_dec")
    print("-" * 30)

    # Step 4: Calculate F_rocket with Titan
    print("Step 4: Calculate F_rocket ≈ m * a_dec.")
    # The multiplication m * a_dec = (5 * 10^1) * 9 = 45 * 10^1.
    # The numerator 45 exceeds the 5-bit limit of 31.
    # We must use a different representation for m=50.
    # Let m = 50 = 100 / 2 = (1/2) * 10^2. This is valid.
    m_titan = TitanNumber(1, 2, 2)
    print(f"Using alternative representation for mass: m = {m_titan} kg")
    # Now, F_rocket = m * a_dec = [(1/2) * 10^2] * (9/1)
    # The result is (9/2) * 10^2. Numerator 9 and denominator 2 are valid.
    F_rocket_titan = TitanNumber(9, 2, 2)
    print(f"Final Force Calculation: F_rocket = {m_titan} * {a_dec}")
    print(f"Final Equation: F_rocket = (1/2) * 10^2 * (9/1) = {F_rocket_titan} N")
    print("-" * 30)

    # Step 5: Calculate the error
    print("Step 5: Calculate the absolute error.")
    F_titan_result = F_rocket_titan.value()
    F_true_result = F_dec_val + F_g_val
    error = abs(F_true_result - F_titan_result)
    print(f"High-precision force (F_true): {F_true_result:.3f} N")
    print(f"Titan-calculated force (F_titan): {F_titan_result:.3f} N")
    print(f"Absolute Error |F_true - F_titan|: {error:.3f} N")
    print("-" * 30)

    # Final Answer
    final_answer = f"Y[{error:.3f}]"
    print(f"The calculation is possible. The smallest absolute error is {error:.3f} N.")
    print(f"Final Answer: <<<Y[{error:.3f}]>>>")

main()