import math

def get_prime_factorization(n):
    """
    Calculates the prime factorization of a positive integer n.
    Returns a dictionary where keys are prime factors and values are their exponents.
    """
    if n <= 0:
        raise ValueError("Input must be a positive integer.")
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while temp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def main():
    """
    Calculates the cardinality of T_l for a given l and prints the steps.
    """
    l = 60
    
    print(f"For l = {l}:")

    # Part A
    U_l_cardinality = l
    print(f"A) |U_l| = l = {U_l_cardinality}")

    # Part B
    if l == 1:
        s = 0
        e = []
        d = 1
        delta_s0 = 1
        P = 1
        T_l_cardinality = 1
        print(f"B) l = 1, so s = 0, d = 1.")
        print(f"   |T_l| = P - 1 + d*delta_s0 = 1 - 1 + 1*1 = 1")

    else:
        factors = get_prime_factorization(l)
        e = list(factors.values())
        s = len(e)
        d = 1 if l % 2 != 0 else 0
        delta_s0 = 0 # since l > 1, s > 0

        # Build the expression string
        p_terms = [f"(2*{exp}+1)" for exp in e]
        p_expr = "*".join(p_terms)
        
        # Calculate P
        P = 1
        for exp in e:
            P *= (2 * exp + 1)
        
        T_l_cardinality = P - 1

        print(f"B) The prime factorization of {l} is {' * '.join([f'{p}^{factors[p]}' for p in factors])}.")
        print(f"   The exponents are e_i = {e}, so s = {s}.")
        print(f"   l is {'even' if d == 0 else 'odd'}, so d = {d}.")
        print(f"   Since s > 0, delta_s0 = 0.")
        print(f"   |T_l| = (product of (2*e_i+1)) - 1 + d*delta_s0")
        print(f"   |T_l| = {p_expr} - 1 + {d}*0")
        p_values = [2 * exp + 1 for exp in e]
        p_values_str = "*".join(map(str, p_values))
        print(f"   |T_l| = {p_values_str} - 1")
        print(f"   |T_l| = {P} - 1 = {T_l_cardinality}")


if __name__ == "__main__":
    main()
