import math

def solve_cardinality(l: int):
    """
    Calculates the cardinalities |U_l| and |T_l| based on the provided definitions
    and prints the results with step-by-step evaluation.

    Args:
        l: A positive integer.
    """
    if not isinstance(l, int) or l < 1:
        print("Error: Input l must be a positive integer.")
        return

    print(f"Solving for l = {l}:")
    print("=" * 25)

    # --- Prime Factorization ---
    factors = {}
    n = l
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    
    exponents = list(factors.values())
    s = len(exponents)

    # --- Part A: |U_l| Calculation ---
    val_A = 2**s
    print("A) Formula: |U_l| = 2^s, where s is the number of distinct prime factors.")
    if l == 1:
        print(f"   l = 1 has s = 0 distinct prime factors.")
    else:
        print(f"   The distinct prime factors of l = {l} are {list(factors.keys())}, so s = {s}.")
    print(f"   Result: |U_{l}| = 2^{s} = {val_A}")
    print("-" * 25)

    # --- Part B: |T_l| Calculation ---
    print("B) Formula for |T_l| is derived from its definition:")
    if l == 1:
        # For l=1, the definition yields |T_1|=1 directly.
        val_B = 1
        print("   For l = 1, the only valid triple is (m,n,λ)=(1,1,1).")
        print(f"   Result: |T_1| = 1")
    else:
        # For l > 1, the formula is (product of (2*e_i + 1)) - 1.
        prod_expr_terms = []
        prod_val_terms = []
        prod_val = 1
        for p, e in factors.items():
            prod_expr_terms.append(f"(2*{e}+1)")
            term_val = 2 * e + 1
            prod_val_terms.append(str(term_val))
            prod_val *= term_val
        
        val_B = prod_val - 1
        
        prod_expr_str = "*".join(prod_expr_terms)
        prod_val_str = "*".join(prod_val_terms)

        print(f"   For l > 1, |T_l| = (Π (2*e_i + 1)) - 1, where e_i are the exponents in the prime factorization.")
        print(f"   The exponents for l = {l} are {exponents}.")
        print(f"   Calculation: |T_{l}| = {prod_expr_str} - 1")
        if len(exponents) > 1:
            print(f"                     = {prod_val_str} - 1")
        print(f"                     = {prod_val} - 1")
        print(f"   Result: |T_{l}| = {val_B}")

# --- Example Execution ---
# You can change the value of `l_to_solve` to test other integers.
l_to_solve = 12
solve_cardinality(l_to_solve)
