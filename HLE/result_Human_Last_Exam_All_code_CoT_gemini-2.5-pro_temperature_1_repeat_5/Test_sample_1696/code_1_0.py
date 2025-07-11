import math
from scipy.special import gamma

def F(z):
    """The Gamma function."""
    return gamma(z)

def G(z):
    """The reciprocal of the Gamma function."""
    if z == 0 or (z < 0 and z == math.floor(z)):
        # Handle poles of Gamma function at 0 and negative integers
        return 0
    return 1 / gamma(z)

def T_val(a, b):
    """Computes the value of the term in the second product."""
    term1 = (2 * a / b) - 1
    term2 = 2**(1 - (2 * a / b))
    term3 = math.sqrt(math.pi)
    term4 = G((2 * a - b) / b)
    return term1 * term2 * term3 * term4

def T_str(a, b):
    """Creates the string representation of the term in the second product."""
    return f"((({2*a}/{b}-1)) * 2^(1-{2*a}/{b}) * sqrt(pi) * G(({2*a-b})/{b}))"

def solve_p_n(n):
    """
    Constructs the formula for P(n), calculates its value, and prints the result.
    """
    if not isinstance(n, int) or n <= 0:
        print("Please provide a positive integer for n.")
        return

    total_value = 1.0
    formula_parts = []

    for b in range(1, n + 1):
        floor_n_b = n // b
        if floor_n_b == 0:
            continue

        b_term_val = 1.0
        b_term_str_parts = []

        if b % 4 != 0:
            # First product condition
            for a in range(1, b):
                if math.gcd(a, b) == 1:
                    b_term_val *= F(a / b)
                    b_term_str_parts.append(f"F({a}/{b})")
        else: # b is divisible by 4
            # Second product condition
            # Note: if gcd(a,b)=1 and 4|b, then a cannot be b/2.
            for a in range(1, b // 2 + 1):
                if math.gcd(a, b) == 1:
                    b_term_val *= T_val(a, b)
                    b_term_str_parts.append(T_str(a, b))

        if b_term_str_parts:
            # Combine parts for this b
            if len(b_term_str_parts) > 1:
                b_term_str = f"({'*'.join(b_term_str_parts)})^({floor_n_b})"
            else:
                b_term_str = f"({b_term_str_parts[0]})^({floor_n_b})"
            
            formula_parts.append(b_term_str)
            total_value *= (b_term_val**floor_n_b)

    final_formula_str = f"P({n}) = {' * '.join(formula_parts)}"
    print(final_formula_str)
    print(f"P({n}) = {total_value}")
    print(f"<<<{total_value}>>>")

# Example calculation for n=4
solve_p_n(4)