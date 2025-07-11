import math

# Helper function for Greatest Common Divisor, needed for fraction simplification
def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    return math.gcd(a, b)

class Fraction:
    """
    A class to represent a fraction in the Titan 5-bit architecture.
    The numerator and denominator must be between 0 and 31.
    """
    def __init__(self, num, den):
        if not (0 <= num <= 31 and 0 <= den <= 31):
            raise ValueError("Numerator and denominator must be 5-bit integers (0-31)")
        if den == 0:
            raise ValueError("Denominator cannot be zero")
        self.num = num
        self.den = den

    def __str__(self):
        return f"{self.num}/{self.den}"

    def __mul__(self, other):
        """
        Multiplies two fractions according to Titan's rules.
        It uses cross-cancellation to simplify before multiplication.
        The resulting numerator and denominator of the simplified fraction must be within the 5-bit range.
        """
        # Simplify before multiplying to keep intermediate values small (cross-cancellation)
        g1 = gcd(self.num, other.den)
        g2 = gcd(other.num, self.den)

        new_num_part1 = self.num // g1
        new_num_part2 = other.num // g2
        
        new_den_part1 = self.den // g2
        new_den_part2 = other.den // g1

        # The final numerator and denominator are calculated
        final_num = new_num_part1 * new_num_part2
        final_den = new_den_part1 * new_den_part2

        if not (0 <= final_num <= 31 and 0 <= final_den <= 31):
            # This would be an error on the Titan computer
            raise ValueError(f"Multiplication result {final_num}/{final_den} is out of 5-bit bounds")

        return Fraction(final_num, final_den)

def solve():
    """
    Solves the Pandora landing problem using the Titan architecture simulation.
    """
    # Step 1 & 2: Calculate true values for g and F (outside Titan)
    G = 6.674e-11
    R_core = 50e3 # m
    R_total = 1000e3 # m
    rho_core = 1200 # kg/m^3
    rho_shell = 300 # kg/m^3
    m_probe = 30 # kg
    h = 500 # m

    # Volume of the core and shell
    V_core = (4/3) * math.pi * R_core**3
    V_shell = (4/3) * math.pi * (R_total**3 - R_core**3)

    # Mass of the core and shell
    M_core = V_core * rho_core
    M_shell = V_shell * rho_shell
    M_pandora = M_core + M_shell

    # Distance from center to probe
    r_probe = R_total + h

    # True gravitational acceleration (g) and force (F)
    g_true = (G * M_pandora) / r_probe**2
    F_true = m_probe * g_true

    # Step 3, 4, 5: Approximate g with a pre-computed 5-bit fraction
    # g_true is approx 0.0838. The fraction 1/12 is 0.0833..., a good approximation.
    g_approx_f = Fraction(1, 12)
    
    # Step 6: Perform the calculation on Titan
    m_probe_f = Fraction(30, 1)
    
    try:
        # This is the single calculation performed by the Titan computer
        F_result_f = m_probe_f * g_approx_f
        
        # Step 7: Calculate error and format output
        F_calculated = F_result_f.num / F_result_f.den
        error = abs(F_calculated - F_true)

        print("Titan Calculation Steps:")
        print(f"1. Approximate Pandora's gravitational acceleration 'g' ({g_true:.4f} m/s^2) with a pre-computed fraction.")
        print(f"   Selected g_frac = {g_approx_f}")
        print("2. Represent the probe's mass 'm' as a fraction.")
        print(f"   Selected m_frac = {m_probe_f}")
        print("3. Perform the multiplication F = m * g on Titan.")
        
        # Final equation as requested
        print("\nFinal Equation:")
        print(f"{F_result_f.num}/{F_result_f.den} = {m_probe_f.num}/{m_probe_f.den} * {g_approx_f.num}/{g_approx_f.den}")
        
        print(f"\nCalculated Force: {F_calculated:.3f} N")
        print(f"True Force: {F_true:.3f} N")
        print(f"Absolute Error: {error:.3f}")
        
        # Final answer format
        final_answer = f"Y[{error:.3f}]"
        print(f"\n<<< {final_answer} >>>")

    except ValueError as e:
        print(f"Calculation failed: {e}")
        print("\n<<< N0 >>>")

solve()