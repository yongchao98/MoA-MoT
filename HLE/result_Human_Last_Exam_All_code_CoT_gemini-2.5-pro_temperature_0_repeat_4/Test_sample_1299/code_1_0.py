import math

def get_prime_factorization(n):
    """
    Returns a dictionary of prime factors and their exponents.
    e.g., for 12, returns {2: 2, 3: 1}
    """
    factors = {}
    d = 2
    temp_n = n
    while d * d <= temp_n:
        while temp_n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp_n //= d
        d += 1
    if temp_n > 1:
        factors[temp_n] = factors.get(temp_n, 0) + 1
    return factors

def solve_cardinality(l):
    """
    Calculates and prints the cardinalities of U_l and T_l.
    """
    if not isinstance(l, int) or l <= 0:
        print("Error: l must be a positive integer.")
        return

    factors = get_prime_factorization(l)
    s = len(factors)
    exponents = list(factors.values())

    # Part A
    # |U_l| = 2^s
    card_U_l = 2**s
    
    # Part B
    # |T_l| = (product of (1+2e_i)) - 1 + delta_{l,1}
    if l == 1:
        card_T_l = 1
        prod_term_str = "1"
        prod_val = 1
    else:
        prod_val = 1
        prod_term_parts = []
        for e in exponents:
            prod_val *= (1 + 2 * e)
            prod_term_parts.append(f"(1+2*{e})")
        prod_term_str = "*".join(prod_term_parts)
        card_T_l = prod_val - 1

    # Print the results in the specified format
    print(f"For l = {l}:")
    print(f"Prime factorization: {' * '.join([f'{p}^{e}' for p, e in factors.items()]) if l > 1 else '1'}")
    print(f"s = {s}, exponents e_i = {exponents}")
    
    # Output for A
    print(f"A) |U_l| = 2^s = 2^{s} = {card_U_l}")
    
    # Output for B
    if l == 1:
        print(f"B) |T_l| = 1")
    else:
        print(f"B) |T_l| = (product of (1+2e_i)) - 1 = {prod_term_str} - 1 = {prod_val} - 1 = {card_T_l}")

# Example usage with a given l, e.g., l = 12
l = 12
solve_cardinality(l)

# Example usage with l = 1
print("\n" + "="*20 + "\n")
solve_cardinality(1)