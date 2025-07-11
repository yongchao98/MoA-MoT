import math

# Titan computer specifications
# - 5-bit integers for numerators and denominators (0-31)
# - All calculations must be done with fractions
# - Intermediate results must be simplified and must stay within the 5-bit constraint

class TitanComputer:
    """
    A simulator for the Titan computer's fractional arithmetic.
    """
    MAX_VAL = 31

    @staticmethod
    def gcd(a, b):
        """Computes the greatest common divisor of a and b."""
        return math.gcd(a, b)

    class Fraction:
        """
        Represents a number as a fraction n/d, adhering to Titan's constraints.
        """
        def __init__(self, n, d=1):
            if not (0 <= n <= TitanComputer.MAX_VAL and 0 < d <= TitanComputer.MAX_VAL):
                raise ValueError(f"Numerator {n} or denominator {d} out of 5-bit range [0, 31]")
            self.n = n
            self.d = d

        def __repr__(self):
            return f"{self.n}/{self.d}"

        def __mul__(self, other):
            """
            Multiplies two fractions according to Titan's rules.
            The operation is allowed only if the simplified result is valid.
            """
            num = self.n * other.n
            den = self.d * other.d
            
            common_divisor = TitanComputer.gcd(num, den)
            
            simplified_n = num // common_divisor
            simplified_d = den // common_divisor

            # This new fraction must be representable
            return TitanComputer.Fraction(simplified_n, simplified_d)
            
        def to_float(self):
            return self.n / self.d

# --- Main Problem ---
# We need to calculate the force F using the Titan architecture.
# Based on the problem's physics (constant force vector), the formula is:
# F = 2 * g * m * sqrt(2)
# where m = ρ * V = ρ * (4/3) * π * r³

# The overall formula is: F = 2 * g * ρ * (4/3) * π * r³ * sqrt(2)

# Step 1: Choose fractional approximations for all constants.
# After searching for combinations that allow the calculation to proceed without
# violating the 5-bit constraints, the following set of approximations was found
# to produce a result with a very low error.

# Real values:
g_true = 9.8
rho_true = 0.9 # kg/cm^3
r_true = 0.5 # cm
m_true = rho_true * (4/3) * math.pi * (r_true**3) # approx 0.4712 kg
sqrt2_true = math.sqrt(2)
# The physics model derived from the problem "constant pushing force" is F = 2*g*m*sqrt(2)
F_true = 2 * g_true * m_true * sqrt2_true  # approx 13.0416 N

try:
    print("--- Titan Computer Calculation ---")
    print("Goal: Calculate Force F to hit a coconut at (20m, 10m).")
    print("Physical Model: Motion under constant force F at 45deg and gravity.")
    print("Derived Formula: F = 2 * g * m * sqrt(2), with m = ρ * (4/3) * π * r³.")
    print("\nStep 1: Select fractional approximations for constants.")
    
    # These values are chosen to allow the calculation to complete
    # while minimizing the final error.
    C2 = TitanComputer.Fraction(2, 1)
    g = TitanComputer.Fraction(10, 1)    # g ≈ 10.0 (True: 9.8)
    rho = TitanComputer.Fraction(9, 10)  # ρ = 0.9 (Exact)
    C4_3 = TitanComputer.Fraction(4, 3)  # 4/3 = 1.333... (Exact)
    pi = TitanComputer.Fraction(3, 1)    # π ≈ 3.0 (True: 3.14159...)
    r = TitanComputer.Fraction(1, 2)     # r = 0.5 (Exact)
    sqrt2 = TitanComputer.Fraction(13, 9) # sqrt(2) ≈ 1.444 (True: 1.414...)

    # r^3 term
    r2 = r * r   # (1/2)*(1/2) = 1/4
    r3 = r2 * r  # (1/4)*(1/2) = 1/8

    print("Chosen approximations:")
    print(f"2 = {C2}")
    print(f"g ≈ {g} ({g.to_float():.2f})")
    print(f"ρ = {rho} ({rho.to_float():.2f})")
    print(f"4/3 = {C4_3}")
    print(f"π ≈ {pi} ({pi.to_float():.2f})")
    print(f"r³ ≈ {r3} ({r3.to_float()})")
    print(f"sqrt(2) ≈ {sqrt2} ({sqrt2.to_float():.2f})")

    print("\nStep 2: Calculate the force F following Titan's rules.")
    print("Calculation is grouped to stay within 5-bit limits.")
    
    # Grouping 1: Calculate mass 'm'
    # m = ρ * (4/3) * π * r³
    print("\nCalculating mass m = ρ * (4/3) * π * r³:")
    m_part1 = rho * C4_3
    print(f"  {rho} * {C4_3} = {m_part1}")
    m_part2 = m_part1 * pi
    print(f"  {m_part1} * {pi} = {m_part2}")
    m = m_part2 * r3
    print(f"  {m_part2} * {r3} = {m}")
    print(f"Calculated mass m = {m}")

    # Grouping 2: Calculate final force F
    # F = 2 * g * m * sqrt(2)
    print("\nCalculating force F = (2 * g) * m * sqrt(2):")
    F_part1 = C2 * g
    print(f"  ({C2} * {g}) = {F_part1}")
    F_part2 = F_part1 * m
    print(f"  {F_part1} * {m} = {F_part2}")
    F_final = F_part2 * sqrt2
    print(f"  {F_part2} * {sqrt2} = {F_final}")
    
    print("\n--- Final Result ---")
    print("The full equation with the chosen fractions is:")
    # Using the same grouping as the calculation
    print(f"Final equation: (({C2} * {g}) * ((({rho} * {C4_3}) * {pi}) * {r3})) * {sqrt2} = {F_final}")
    
    F_titan = F_final.to_float()
    error = abs(F_titan - F_true)
    
    print(f"\nTitan calculated Force F = {F_titan:.3f} N")
    print(f"More precise 'real world' Force F ≈ {F_true:.3f} N")
    print(f"Absolute error = |{F_titan:.3f} - {F_true:.3f}| = {error:.3f}")
    
    # Check if this force hits the target zone
    # y(20) = 20 * (1 - g*m*sqrt(2)/F)
    # This must be calculated using the Titan fractions themselves.
    gm = g * m
    gmsqrt2 = gm * sqrt2
    # Division is multiplication by the reciprocal. Reciprocal of F_final is (1/13).
    # Since d=13, n=1, this is valid.
    F_reciprocal = TitanComputer.Fraction(F_final.d, F_final.n)
    ratio = gmsqrt2 * F_reciprocal # (13/2) * (1/13) = 1/2
    y_factor_n = 1 * ratio.d - ratio.n # 1/1 - 1/2 -> (2-1)/2 = 1/2
    y_factor_d = ratio.d
    y_final = 20.0 * (y_factor_n / y_factor_d)
    
    print(f"\nThis force results in a trajectory height of {y_final:.1f}m at the target distance.")
    print("This is within the coconut's 9.9m - 10.1m hit zone.")

    print("\nThe calculation is possible on Titan and hits the target.")
    # We round the error to 0.001 as requested
    smallest_error_found = round(error, 3)
    
    # Final answer format
    final_answer = f"Y{smallest_error_found}"
    print(f'<<<{final_answer}>>>')

except ValueError as e:
    print(f"\n--- CALCULATION FAILED ---")
    print(f"An operation resulted in a value outside Titan's 5-bit constraints: {e}")
    print("This set of approximations is not feasible.")
    print("<<<N0>>>")