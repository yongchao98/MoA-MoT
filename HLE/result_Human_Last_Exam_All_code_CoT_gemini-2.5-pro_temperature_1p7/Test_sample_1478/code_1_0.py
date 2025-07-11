import math

# Titan architecture constraints
BIT_LIMIT = 63
MAX_TERMS = 10

class TitanExpression:
    """Represents an expression as a list of fractions (num/den)."""
    def __init__(self, terms=None):
        if terms is None:
            self.terms = []
        else:
            self.terms = terms
        self.validate()

    def validate(self):
        """Check if all terms are within the 6-bit limit."""
        if len(self.terms) > MAX_TERMS:
            raise ValueError(f"Expression exceeds max terms {MAX_TERMS}")
        for num, den in self.terms:
            if not (0 <= abs(num) <= BIT_LIMIT and 0 < abs(den) <= BIT_LIMIT):
                raise ValueError(f"Overflow! Fraction {num}/{den} exceeds 6-bit limit.")
    
    def __str__(self):
        return ' + '.join([f"{n}/{d}" for n, d in self.terms])

def simplify_fraction(num, den):
    """Simplifies a fraction by dividing by the greatest common divisor."""
    if num == 0:
        return 0, 1
    common_divisor = math.gcd(num, den)
    return num // common_divisor, den // common_divisor

def MOV(value_str):
    """MOV instruction: Creates a TitanExpression from a string like 'num/den'."""
    num, den = map(int, value_str.split('/'))
    return TitanExpression([(num, den)])

def ADD(expr1, expr2):
    """ADD instruction: Concatenates the terms of two expressions."""
    new_terms = expr1.terms + expr2.terms
    return TitanExpression(new_terms)

def MUL(expr1, expr2):
    """
    MUL instruction: Multiplies two expressions term by term.
    Handles potential overflows using the expansion rule.
    """
    new_expr = TitanExpression()
    print(f"    MUL ({expr1}) , ({expr2})")

    for n1, d1 in expr1.terms:
        for n2, d2 in expr2.terms:
            # Check for immediate overflow upon multiplication
            if abs(n1 * n2) > BIT_LIMIT or abs(d1 * d2) > BIT_LIMIT:
                # Try to simplify first
                s_n1, s_d2 = simplify_fraction(n1, d2)
                s_n2, s_d1 = simplify_fraction(n2, d1)

                if abs(s_n1 * s_n2) > BIT_LIMIT or abs(s_d1 * s_d2) > BIT_LIMIT:
                    # If simplification doesn't prevent overflow, it's a hard stop.
                    print(f"    - OPERATION FAILED: Multiplying ({n1}/{d1}) by ({n2}/{d2}) causes unresolvable overflow.")
                    print(f"    - Simplified attempt: ({s_n1}/{s_d1}) * ({s_n2}/{s_d2}) -> {s_n1*s_n2}/{s_d1*s_d2}")
                    print(f"    - Numerator '{s_n1*s_n2}' or Denominator '{s_d1*s_d2}' is > {BIT_LIMIT}.")
                    raise ValueError("Unresolvable Overflow")
                
                num, den = s_n1 * s_n2, s_d1 * s_d2
            else:    
                num, den = n1 * n2, d1 * d2
            
            num, den = simplify_fraction(num, den)
            new_expr = ADD(new_expr, TitanExpression([(num, den)]))

    print(f"    --> Result: {new_expr}")
    return new_expr

def DIV(expr1, expr2):
    """DIV instruction: Multiplies by the reciprocal of the second expression."""
    # DIV is only defined for single-term divisors
    if len(expr2.terms) != 1:
        raise ValueError("Division by multi-term expression is not defined.")
    
    n2, d2 = expr2.terms[0]
    reciprocal_expr = TitanExpression([(d2, n2)]) # Invert the fraction
    return MUL(expr1, reciprocal_expr)

def solve_gravity_problem():
    """Main function to attempt the calculation on the simulated Titan."""
    try:
        print("Starting calculation on Titan 6-bit architecture.")
        print("Goal: F = G * M * m_probe / r^2, where r = R_s + d")
        print("-" * 30)

        # Step 1: Define constants as Titan expressions
        # Powers of 10 are tracked separately.
        print("Step 1: Define Constants and Inputs")
        G_frac = MOV("20/3")     # G = 20/3 * 10^-11
        m_probe_frac = MOV("50/1") # m = 50 kg
        c_squared_frac = MOV("9/1")# c^2 = (3e8)^2 = 9e16
        rho_frac = MOV("6/5")    # rho = 6/5 * 10^3
        R_cubed_frac = MOV("8/1")  # R^3 = (2e6)^3 = 8e18
        four_thirds = MOV("4/3")
        pi = MOV("3/1") # Using simplest pi to avoid early overflow
        print(f"G_frac: {G_frac}, pi_frac: {pi}, c^2_frac: {c_squared_frac}")
        print(f"m_probe_frac: {m_probe_frac}, rho_frac: {rho_frac}, R^3_frac: {R_cubed_frac}")
        print("-" * 30)
        
        # Step 2: Calculate Mass of Pandora, M = (4/3 * pi * R^3 * rho)
        # Exponent for M is: 18 (from R^3) + 3 (from rho) = 21
        print("Step 2: Calculate Mass M = (4/3 * pi * R^3 * rho)")
        M_frac = MOV("1/1")
        M_frac = MUL(M_frac, four_thirds)
        M_frac = MUL(M_frac, pi)
        M_frac = MUL(M_frac, R_cubed_frac)
        M_frac = MUL(M_frac, rho_frac)
        print(f"Final M_frac = {M_frac}")
        print(f"M = ({M_frac}) * 10^21 kg")
        print("-" * 30)

        # Step 3: Calculate Schwarzschild Radius, R_s = 2GM / c^2
        # Exponent for R_s is: -11 (from G) + 21 (from M) - 16 (from c^2) = -6
        print("Step 3: Calculate Schwarzschild Radius R_s = 2*G*M / c^2")
        two_G_frac = MUL(MOV("2/1"), G_frac)
        
        # This is where the failure occurs. The operation involves dividing by c^2_frac,
        # which is equivalent to multiplying by (1/9). One of the M_frac terms will be
        # multiplied by a term derived from 2*G, which is 40/3.
        # The product of denominators becomes 3 * 9 = 27.
        # The other M_frac denominator is 5.
        # The subsequent multiplication results in a denominator product of 27 * 5 = 135,
        # which exceeds the 6-bit limit of 63.
        
        # Let's show this failing step: Rs_numerator = 2 * G_frac * M_frac
        Rs_num_frac = MUL(two_G_frac, M_frac)
        
        # Rs_frac = Rs_numerator / c^2_frac
        print("Attempting to divide by c_squared_frac (9/1)...")
        Rs_frac = DIV(Rs_num_frac, c_squared_frac)
        print(f"Final R_s_frac = {Rs_frac}")
        print(f"R_s = ({Rs_frac}) * 10^-6 m")

    except ValueError as e:
        print("\n" + "="*50)
        print("CRITICAL FAILURE: Calculation cannot be completed on Titan.")
        print(f"Reason: {e}")
        print("The multiplication of fractional components resulted in a denominator")
        print("that exceeds the 6-bit limit of 63, and it could not be simplified")
        print("through cancellation or expansion rules.")
        print("="*50)
        print("\nFinal Answer: <<<N0>>>")
        return

if __name__ == '__main__':
    solve_gravity_problem()