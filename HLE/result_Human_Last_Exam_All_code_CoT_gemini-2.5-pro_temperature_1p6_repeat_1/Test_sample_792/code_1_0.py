import math

def main():
    """
    Calculates the force required for the monkey to hit the lion
    using the computational rules of the Titan architecture.
    """

    # Helper function for Greatest Common Divisor
    def gcd(a, b):
        while b:
            a, b = b, a % b
        return a

    # A class to represent and print fractions for clarity
    class Fraction:
        def __init__(self, num, den=1):
            if den == 0:
                raise ZeroDivisionError
            self.num = num
            self.den = den

        def __repr__(self):
            return f"{self.num}/{self.den}"

    def multiply(f1, f2):
        """
        Simulates Titan's 5-bit fractional multiplication.
        It simplifies fractions before multiplying to avoid overflow.
        Returns a new Fraction object or raises an error if the operation fails.
        """
        num1, den1 = f1.num, f1.den
        num2, den2 = f2.num, f2.den

        # Check if inputs are valid 5-bit integers
        for val in [num1, den1, num2, den2]:
            if not (0 <= val <= 31):
                raise ValueError(f"Input value {val} is not a 5-bit integer.")
        
        # Cross-simplify fractions
        common_divisor1 = gcd(num1, den2)
        num1_s = num1 // common_divisor1
        den2_s = den2 // common_divisor1

        common_divisor2 = gcd(num2, den1)
        num2_s = num2 // common_divisor2
        den1_s = den1 // common_divisor2

        # Perform multiplication on simplified parts
        res_num = num1_s * num2_s
        res_den = den1_s * den2_s

        # Check if the result is within 5-bit limits
        if res_num > 31 or res_den > 31:
            raise ValueError(f"Calculation failed: ({num1}/{den1}) * ({num2}/{den2}) -> {res_num}/{res_den}. Values exceed 31.")
        
        return Fraction(res_num, res_den)

    # --- Problem Setup ---
    # Define initial parameters as Titan fractions
    print("Step 1: Define initial physical parameters as Titan fractions.")
    r = Fraction(1, 2)    # radius = 0.5 cm
    rho = Fraction(9, 10) # density = 0.9 kg/cm^3
    d = Fraction(20, 1)   # distance = 20 m
    h = Fraction(10, 1)   # height = 10 m
    
    # Approximations for constants that allow calculation within constraints
    print("Step 2: Select workable approximations for g and π.")
    g = Fraction(10, 1)   # g ≈ 10 m/s^2
    pi = Fraction(28, 9)  # π ≈ 3.111...
    const_4_3 = Fraction(4, 3)
    
    print(f"Using: r={r}, ρ={rho}, d={d}, h={h}, g={g}, π={pi}\n")

    # --- Calculation ---
    print("Step 3: Calculate the mass (m) of the rock.")
    print("m = (4/3) * π * r³ * ρ")
    # r² = r * r
    r2 = multiply(r, r)
    print(f"r² = {r} * {r} = {r2}")
    # r³ = r² * r
    r3 = multiply(r2, r)
    print(f"r³ = {r2} * {r} = {r3}")
    # Partial mass calculation to manage complexity
    m_part1 = multiply(const_4_3, rho)
    print(f"(4/3) * ρ = {const_4_3} * {rho} = {m_part1}")
    m_part2 = multiply(m_part1, r3)
    print(f"({m_part1}) * r³ = {m_part1} * {r3} = {m_part2}")
    m = multiply(m_part2, pi)
    print(f"m = ({m_part2}) * π = {m_part2} * {pi} = {m}\n")
    
    print("Step 4: Calculate the required force (F).")
    print("F = (d/h) * g * m")
    # d_over_h = d / h = d * (1/h)
    h_inv = Fraction(h.den, h.num)
    d_over_h = multiply(d, h_inv)
    print(f"d/h = {d} / {h} = {d_over_h}")
    # Force calculation
    F_part1 = multiply(d_over_h, g)
    print(f"(d/h) * g = {d_over_h} * {g} = {F_part1}")
    F_final = multiply(F_part1, m)
    
    # --- Final Result ---
    print("\n--- Final Calculation Summary ---")
    final_equation = f"F = ({d.num}/{d.den} / {h.num}/{h.den}) * ({g.num}/{g.den}) * ({const_4_3.num}/{const_4_3.den}) * ({pi.num}/{pi.den}) * ({r.num}/{r.den})³ * ({rho.num}/{rho.den})"
    print(final_equation)
    print(f"The calculated force is {F_final.num}/{F_final.den} Newtons.")

    # Calculate the error
    F_calculated = F_final.num / F_final.den
    F_true = (20 * 9.8 / 10) * (4/3 * math.pi * (0.5**3) * 0.9)
    error = abs(F_calculated - F_true)
    print(f"\nThis corresponds to a decimal value of {F_calculated:.4f} N.")
    print(f"The 'true' physical value is ~{F_true:.4f} N.")
    print(f"The absolute error is {error:.4f} N.")


if __name__ == '__main__':
    main()

# Calculate the error for the final answer tag
F_calculated_val = 28/3
F_true_val = 2.94 * math.pi
error_val = abs(F_calculated_val - F_true_val)
# print(f"<<<Y[{error_val:.6f}]>>>")
# The above line is commented out but used to generate the final answer tag.
# It's here for transparency on how the error value was obtained.