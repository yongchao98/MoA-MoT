import math

# Titan architecture constraints
MAX_VAL = 63
MAX_TERMS = 10

class Fraction:
    """Represents a single fraction with a numerator and denominator."""
    def __init__(self, n, d=1):
        if not (0 <= n <= MAX_VAL and 1 <= d <= MAX_VAL):
            # This check is for constants, simplification handles large values.
            raise ValueError(f"Numerator/Denominator {n}/{d} exceeds 6-bit limit.")
        self.n = n
        self.d = d

    def __repr__(self):
        return f"{self.n}/{self.d}"

class Expression:
    """Represents an expression in a Titan register."""
    def __init__(self, terms, exp=0):
        if len(terms) > MAX_TERMS:
            raise ValueError(f"Expression exceeds the {MAX_TERMS} term limit.")
        self.terms = terms
        self.exp = exp

    def __repr__(self):
        term_str = " + ".join(map(str, self.terms))
        return f"({term_str}) * 10^{self.exp}"

def simplify_integer(n):
    """Simplifies an integer > MAX_VAL into a sum of valid fractions."""
    if n <= MAX_VAL:
        return [Fraction(n)]
    else:
        # Recursively break down the number
        # e.g., 150 -> [63/1] + simplify_integer(87) -> [63/1, 63/1, 24/1]
        return [Fraction(MAX_VAL)] + simplify_integer(n - MAX_VAL)

def simplify_fraction(n, d):
    """Simplifies a fraction N/D into a Titan Expression."""
    common = math.gcd(n, d)
    n //= common
    d //= common

    if n <= MAX_VAL and d <= MAX_VAL:
        return [Fraction(n, d)]

    whole_part = n // d
    remainder = n % d

    terms = simplify_integer(whole_part)
    if remainder > 0:
        # Remainder itself can't be > MAX_VAL if d <= MAX_VAL
        terms.append(Fraction(remainder, d))
    return terms

def multiply_expressions(expr1, expr2):
    """Multiplies two expressions, term by term."""
    new_terms = []
    new_exp = expr1.exp + expr2.exp
    
    print(f"Multiplying:\n  expr1 = {expr1}\n  expr2 = {expr2}")

    for t1 in expr1.terms:
        for t2 in expr2.terms:
            n = t1.n * t2.n
            d = t1.d * t2.d
            simplified_product = simplify_fraction(n, d)
            new_terms.extend(simplified_product)
            if len(new_terms) > MAX_TERMS:
                print(f"--> FAILURE: Term count is {len(new_terms)}, exceeding the limit of {MAX_TERMS}.")
                return None

    print(f"--> SUCCESS: Result has {len(new_terms)} terms.")
    return Expression(new_terms, new_exp)

def add_expressions(expr1, expr2):
    """Adds two expressions. Assumes common exponent for simplicity."""
    if expr1.exp != expr2.exp:
        raise ValueError("Cannot add expressions with different exponents in this simplified simulation.")
    
    # This part is complex (finding common denominators) and not needed to prove the point.
    # The main issue is the multiplicative explosion of terms.
    total_terms = expr1.terms + expr2.terms
    if len(total_terms) > MAX_TERMS:
        print(f"--> FAILURE: Term count for addition is {len(total_terms)}, exceeding limit.")
        return None
    return Expression(total_terms, expr1.exp)

def run_calculation():
    """Simulates the gravity calculation on Titan."""
    print("--- Titan Feasibility Study ---")
    print("Goal: Calculate gravitational force on a probe.")
    print(f"Architecture limits: Numerator/Denominator <= {MAX_VAL}, Terms per register <= {MAX_TERMS}\n")

    # Step 1: Define constants as Titan Expressions
    print("Step 1: Defining initial constants...")
    try:
        G = Expression([Fraction(20, 3)], exp=-11)  # G ~ 6.67e-11
        m = Expression([Fraction(50, 1)])          # probe mass = 50 kg
        
        # We need G*m for the force calculation.
        print("\nStep 2: Calculate G * m")
        Gm_expr = multiply_expressions(G, m)
        if Gm_expr is None: return "N0"
        
        # Result of 20/3 * 50/1 is 1000/3.
        # simplify_fraction(1000, 3): 1000//3=333, 1000%3=1 -> simplify_integer(333) + 1/3
        # simplify_integer(333): 63+63+63+63+63+18.
        # Total terms for G*m: 6 (from 333) + 1 (from 1/3) = 7 terms.
        print(f"Resulting expression for G*m: {Gm_expr}\n")

        # Step 3: Calculate Pandora's Mass (M)
        print("Step 3: Calculate Pandora's Mass M = rho * (4/3) * pi * R^3")
        rho = Expression([Fraction(12, 1)], exp=2)   # 1200 kg/m^3
        four_thirds = Expression([Fraction(4, 3)])
        pi = Expression([Fraction(22, 7)])
        R_cubed = Expression([Fraction(8, 1)], exp=18) # (2e6 m)^3

        # M_part1 = rho * 4/3 = (12e2) * (4/3) = 16e2
        M_part1 = multiply_expressions(rho, four_thirds)
        if M_part1 is None: return "N0"
        print(f"Intermediate M (rho * 4/3): {M_part1}")

        # M_part2 = M_part1 * pi = (16e2) * (22/7) = 352/7 * e2
        # simplify_fraction(352,7): 352//7=50, 352%7=2 -> [50/1, 2/7] (2 terms)
        M_part2 = multiply_expressions(M_part1, pi)
        if M_part2 is None: return "N0"
        print(f"Intermediate M (... * pi): {M_part2}")

        # M = M_part2 * R^3
        M_expr = multiply_expressions(M_part2, R_cubed)
        if M_expr is None: return "N0"

        # The multiplication is (50/1 + 2/7) * (8/1) -> 400/1 + 16/7
        # simplify_fraction(400,1) -> simplify_integer(400) -> 63+63+63+63+63+63+22 (7 terms)
        # simplify_fraction(16,7) -> 16/7 (1 term)
        # Total terms for M: 7 + 1 = 8 terms.
        print(f"Resulting expression for M: {M_expr}\n")

        # Step 4: Attempt to calculate the numerator of the force equation F ~ (G*m) * M / d^2
        # We only need to check if we can compute the numerator, (G*m) * M.
        print("Step 4: Attempt to multiply (G*m) expression by M expression.")
        # Gm_expr has 7 terms. M_expr has 8 terms.
        # The result would have 7 * 8 = 56 sub-results to be simplified and summed.
        # Each sub-result might be more than one term.
        # The total number of terms will vastly exceed the 10-term limit.
        print(f"Expression for (G*m) has {len(Gm_expr.terms)} terms.")
        print(f"Expression for M has {len(M_expr.terms)} terms.")
        print("Multiplying them would result in at least 7 * 8 = 56 terms.")
        print("This exceeds the register's 10-term capacity.")
        
        force_numerator = multiply_expressions(Gm_expr, M_expr)
        if force_numerator is None:
            print("\nConclusion: The calculation is not feasible with the Titan architecture.")
            return "N0"

    except ValueError as e:
        print(f"An error occurred: {e}")
        return "N0"

    # This part will not be reached
    return "Y[...]"

# Execute the simulation and print the final verdict.
final_answer = run_calculation()
print(f"\nFinal Answer: {final_answer}")
print(f"\n<<<{final_answer}>>>")
