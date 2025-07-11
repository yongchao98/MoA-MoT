import math

# Titan computer constants and rules
MAX_VAL = 31

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return a

class TitanFraction:
    """A class to represent and operate on Titan 5-bit fractions."""
    def __init__(self, n, d=1, name=""):
        if not (0 <= n <= MAX_VAL and 1 <= d <= MAX_VAL):
            raise ValueError(f"Invalid 5-bit fraction: {n}/{d}")
        self.n = n
        self.d = d
        self.name = name

    def __repr__(self):
        return f"{self.n}/{self.d}"

    def value(self):
        return self.n / self.d

def multiply(f1, f2, step_name):
    """
    Multiplies two TitanFractions, respecting the 5-bit constraint.
    Returns a new TitanFraction or None if the operation is invalid.
    """
    n1, d1 = f1.n, f1.d
    n2, d2 = f2.n, f2.d
    
    # Simplify before multiplying by finding common divisors
    g1 = gcd(n1, d2)
    n1_s, d2_s = n1 // g1, d2 // g1
    
    g2 = gcd(n2, d1)
    n2_s, d1_s = n2 // g2, d1 // g2
    
    res_n = n1_s * n2_s
    res_d = d1_s * d2_s
    
    print(f"Step: {step_name}")
    print(f"  Multiplying {f1.name} ({f1}) by {f2.name} ({f2})")
    
    if res_n > MAX_VAL or res_d > MAX_VAL:
        print(f"  Result ({res_n}/{res_d}) is INVALID (numerator or denominator > {MAX_VAL}).")
        return None
    
    result_frac = TitanFraction(res_n, res_d)
    print(f"  Result is {result_frac}. This is a valid intermediate fraction.")
    return result_frac

def solve():
    """
    Main function to perform the calculation based on the plan.
    """
    print("Starting Titan calculation for gravitational force on Pioneer.\n")
    print("### Step 1: Define constants as 5-bit fractions ###")
    
    # Using approximations that enable calculation. The key is to find a set of
    # fractions whose factors cancel out nicely.
    # G ≈ 6.67e-11. We choose 21/3 = 7 to enable cancellation with pi.
    # pi ≈ 3.14. We choose 22/7 to enable cancellation with G.
    G = TitanFraction(21, 3, name="G")
    FOUR_THIRDS = TitanFraction(4, 3, name="4/3")
    PI = TitanFraction(22, 7, name="π")
    RHO_SHELL = TitanFraction(3, 1, name="ρ_shell_num") # For 300 kg/m^3
    R_TOTAL = TitanFraction(1, 1, name="R_total_num") # For 1,000,000 m
    M_PROBE = TitanFraction(30, 1, name="m_probe")
    
    # We will compute the numerical part of F ≈ G * (4/3) * π * ρ_shell * R_total * m_probe
    # Powers of 10: G(10^-11), ρ(10^2), R(10^6) -> Total: 10^-3
    
    print(f"Approximations chosen:")
    print(f"  G         ≈ {G.value():.2f}  (using {G})")
    print(f"  π         ≈ {PI.value():.2f}  (using {PI})")
    print(f"  4/3       = {FOUR_THIRDS.value():.2f} (using {FOUR_THIRDS})")
    print(f"  ρ_shell   = 300   (using {RHO_SHELL} x 10^2)")
    print(f"  R_total   = 1e6   (using {R_TOTAL} x 10^6)")
    print(f"  m_probe   = 30    (using {M_PROBE})\n")

    print("### Step 2: Perform multiplication step-by-step ###")
    
    # Order of operations is crucial. We group terms that cancel.
    # F_num = [G * (4/3) * PI * RHO_SHELL] * R_TOTAL * M_PROBE
    
    # Let's combine G and PI first, as they are chosen to cancel.
    # (21/3) * (22/7) -> (7/1) * (22/7) -> 22/1
    temp_G = TitanFraction(7, 1, "G_simplified") # From 21/3
    res1 = multiply(temp_G, PI, "G * PI")
    if res1 is None: return "N0"
    res1.name = "(G*π)"

    # Next, multiply by 4/3 and rho_shell, which also cancel.
    # (4/3) * (3/1) -> 4/1
    res2 = multiply(FOUR_THIRDS, RHO_SHELL, "(4/3) * ρ_shell")
    if res2 is None: return "N0"
    res2.name = "((4/3)*ρ)"

    # Now we have (22/1) * (4/1)
    res3 = multiply(res1, res2, "(G*π) * ((4/3)*ρ)")
    if res3 is None:
        # This is the expected failure point for most approximations.
        # 22 * 4 = 88, which is > 31.
        print("\nCalculation failed. It is not possible to compute the force with these constraints.")
        return "N0"
        
    res3.name = "(G*π*(4/3)*ρ)"

    # If it had succeeded, we would continue...
    res4 = multiply(res3, R_TOTAL, "Multiply by R_total")
    if res4 is None: return "N0"
    res4.name = "(...*R_total)"
    
    final_res = multiply(res4, M_PROBE, "Multiply by m_probe")
    if final_res is None: return "N0"

    # This part will not be reached due to the failure.
    # final_value = final_res.value() * 1e-3 # Apply power of 10
    # print(f"\nFinal numerical part: {final_res}")
    # print(f"Final force F = {final_res.value()} * 10^-3 = {final_value:.3f} N")
    return "Y"


final_answer = solve()
print(f"\nFinal Answer: {final_answer}")
# The problem asks for the error if the answer is Yes. Since the calculation
# fails, the format is just N0.
print("\nSince the calculation failed at an intermediate step due to the 5-bit constraint, the final answer is N0.")
print("The equation to be solved is F = G * (4/3) * pi * rho_shell * R_total * m_probe")
print("Using fractions: F_num = (21/3) * (4/3) * (22/7) * (3/1) * (1/1) * (30/1)")
print("A possible calculation path is: ((21/3) * (22/7)) * ((4/3) * (3/1)) * (30/1)")
print("First term: (21/3) * (22/7) = (7/1) * (22/7) = 22/1")
print("Second term: (4/3) * (3/1) = 4/1")
print("Intermediate product: (22/1) * (4/1) = 88/1")
print("The numerator 88 exceeds the 5-bit limit of 31, so the calculation fails.")
