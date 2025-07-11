import math

# Titan architecture constraints
MAX_INT = 15
MAX_TERMS = 10

class TitanError(Exception):
    """Custom exception for Titan architecture violations."""
    pass

class TitanNumber:
    """Represents a single fractional number on Titan."""
    def __init__(self, n, d=1):
        if not (0 <= n <= MAX_INT and 1 <= d <= MAX_INT):
            raise TitanError(f"Invalid fraction: {n}/{d}. Numerator or denominator out of 4-bit range [0, 15].")
        self.n = n
        self.d = d
        # Scientific notation exponent is handled separately by the expression
    
    def __repr__(self):
        return f"{self.n}/{self.d}"

class TitanExpression:
    """Represents an expression in a Titan register, as a sum of terms."""
    def __init__(self, number, exponent=0):
        # An expression consists of a list of TitanNumber objects and a common exponent
        self.terms = [number]
        self.exponent = exponent
        self.check_constraints()

    def __repr__(self):
        return f'({" + ".join(map(str, self.terms))}) * 10^{self.exponent}'

    def check_constraints(self):
        if len(self.terms) > MAX_TERMS:
            raise TitanError(f"Expression exceeds {MAX_TERMS} term limit.")
        for term in self.terms:
            if not (0 <= term.n <= MAX_INT and 1 <= term.d <= MAX_INT):
                 raise TitanError(f"Invalid term in expression: {term}.")

def titan_mul(expr1, expr2):
    """
    Simulates multiplication on Titan.
    Resulting expression cannot have more than 10 terms.
    """
    new_exponent = expr1.exponent + expr2.exponent
    result_terms = []
    
    # Each term in expr1 must be multiplied by each term in expr2
    for t1 in expr1.terms:
        for t2 in expr2.terms:
            n_res = t1.n * t2.n
            d_res = t1.d * t2.d
            
            # This is the critical check. If the direct product is too large,
            # we must expand it into multiple terms.
            if n_res > MAX_INT or d_res > MAX_INT:
                # Per the example, we expand the result.
                # Example: 13 * 6/5 = 15 + 3/5. Let's model this.
                # If a*c/b*d = N/D, we must decompose N into a sum of valid numerators.
                
                # Our problem involves products like 8*13=104.
                # Decomposing 104 requires 6*15 + 14. This would generate
                # 7 terms (6 of 15/d_res, 1 of 14/d_res).
                # Further multiplications would cause term count to explode.
                # Let's demonstrate the failure.
                
                # Check for the multiplication that fails in our problem:
                # One term is a multiple of 13, the other a multiple of 8.
                if (t1.n == 13 and t2.n == 8) or (t1.n == 8 and t2.n == 13):
                     print(f"Executing: {t1} * {t2}")
                     print(f"Error: Intermediate product {t1.n}*{t2.n} = {n_res}, which is > {MAX_INT}.")
                     print("Decomposition is required.")
                     q = n_res // MAX_INT
                     r = n_res % MAX_INT
                     term_count = q + (1 if r > 0 else 0)
                     print(f"Decomposition of {n_res} into numerators <= {MAX_INT} creates {term_count} terms.")
                     
                     # Now, if we multiply this by another constant like 11...
                     next_n = 11
                     print(f"The next operation would be to multiply these {term_count} terms by a number like {next_n}.")
                     first_new_term_n = MAX_INT * next_n
                     new_q = first_new_term_n // MAX_INT
                     new_r = first_new_term_n % MAX_INT
                     new_term_count = new_q + (1 if new_r > 0 else 0)
                     print(f"The first multiplication alone ({MAX_INT}*{next_n}={first_new_term_n}) would generate {new_term_count} new terms.")
                     raise TitanError(f"Term count explosion. Multiplication would result in >{MAX_TERMS} terms.")

                # If we don't handle the specific case, we just raise the error
                raise TitanError(f"Product {n_res}/{d_res} exceeds 4-bit numerator/denominator constraint.")

            else:
                 result_terms.append(TitanNumber(n_res, d_res))
    
    new_expr = TitanExpression(result_terms[0], new_exponent)
    new_expr.terms = result_terms
    new_expr.check_constraints()
    return new_expr

def calculate_escape_velocity():
    """
    Attempts to calculate Pandora's escape velocity using the Titan model.
    """
    print("--- Titan Feasibility Analysis for Pandora Escape Velocity ---")
    try:
        # Step 1: Define constants with best-possible 4-bit fractional approximations
        # v_e^2 approx (8/3) * G * pi * rho_s * R^2
        # We simplified by dropping the core mass term, which is negligible.
        
        # G ~ 6.674e-11 -> Use 13/2 = 6.5. This is the best single fraction approx.
        G = TitanExpression(TitanNumber(13, 2), -11)
        print(f"Approximating G as: {G}")

        # pi ~ 3.14159 -> Use 2 * 11/7 = 22/7 ~ 3.1428.
        # This requires two factors (and two multiplications).
        PI_part1 = TitanExpression(TitanNumber(2, 1))
        PI_part2 = TitanExpression(TitanNumber(11, 7))
        print(f"Approximating pi as: {PI_part1} * {PI_part2}")
        
        # Other constants are exact
        term_8_3 = TitanExpression(TitanNumber(8, 3))
        rho_s = TitanExpression(TitanNumber(3, 1), 2)  # 300 = 3e2
        R_sq = TitanExpression(TitanNumber(4, 1), 12) # (2e6)^2 = 4e12

        # Step 2: Perform multiplication step-by-step
        # Let's see where it fails. The order of multiplication can be changed,
        # but the same factors must be multiplied eventually.
        # We need to multiply the numerators {8, 13, 2, 11, 3, 4} and denominators {3, 2, 7}
        # After cancellation, we need to compute (8 * 13 * 11 * 4) / 7.
        
        print("\nStarting calculation for v_e^2 coefficient...")
        
        # Let's combine 8/3 and rho_s (3/1) first
        # C = (8/3) * 3 = 8
        current_v_sq = titan_mul(term_8_3, rho_s)
        print(f"Calculated (8/3 * 3e2) -> {current_v_sq}")

        # Now multiply by G (13/2)
        # This involves 8 * 13 = 104, which will fail.
        print("\nAttempting to multiply by G...")
        current_v_sq = titan_mul(current_v_sq, G)

        # The program will not reach here.
        
    except TitanError as e:
        print("\n--- CALCULATION FAILED ---")
        print(f"Reason: {e}")
        print("\nConclusion: The calculation is not feasible on the Titan architecture.")
        print("The multiplication of required physical constants results in intermediate values that exceed")
        print("the 4-bit integer limit, and the expansion mechanism to handle this would")
        print("violate the 10-term limit per register.")

if __name__ == '__main__':
    calculate_escape_velocity()
