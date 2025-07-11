import math

def get_prime_factorization_exponents(n):
    """
    Returns a list of exponents in the prime factorization of n.
    For example, for n=12=2^2 * 3^1, this returns [2, 1].
    """
    if n <= 0:
        return []
    exponents = []
    
    count = 0
    while n % 2 == 0:
        count += 1
        n //= 2
    if count > 0:
        exponents.append(count)

    d = 3
    while d * d <= n:
        count = 0
        while n % d == 0:
            count += 1
            n //= d
        if count > 0:
            exponents.append(count)
        d += 2
    
    if n > 1:
        exponents.append(1)
        
    return exponents

def solve_for_l(l):
    """
    Calculates and prints the cardinalities of T_l and U_l for a given integer l.
    It breaks down the calculation as requested.
    """
    print(f"--- For l = {l} ---")

    if l == 1:
        T_val = 1
        # d=1 since l=1 is odd
        d = 1
        U_val = T_val + d
        print(f"A)[For l=1, |U_1| = |T_1| + d = 1 + 1 = {U_val}]")
        print(f"B)[For l=1, the only solution is (m,n,lambda)=(1,1,1), so |T_1| = {T_val}]")
        return

    # For l > 1
    exponents = get_prime_factorization_exponents(l)
    
    # Calculate |T_l| = tau(l^2) - 1
    tau_l_sq_terms = [2 * e + 1 for e in exponents]
    
    # Calculate tau(l^2) by multiplying the terms
    tau_l_sq = 1
    for term in tau_l_sq_terms:
        tau_l_sq *= term
        
    T_val = tau_l_sq - 1

    # Calculate |U_l| = |T_l| + d
    d = 0 if l % 2 == 0 else 1
    U_val = T_val + d
    
    # Format the expression strings for B
    b_expr_parts = [f"(2*{e}+1)" for e in exponents]
    b_expr_str = " * ".join(b_expr_parts)
    b_terms_eval_str = " * ".join(map(str, tau_l_sq_terms))

    print(f"B)[|T_{l}| = ( {b_expr_str} ) - 1 = ( {b_terms_eval_str} ) - 1 = {T_val}]")
    
    # Format the expression strings for A
    if d == 0:  # l is even and > 1
        print(f"A)[|U_{l}| = |T_{l}| + 0 = {T_val}]")
    else:  # l is odd and > 1
        print(f"A)[|U_{l}| = |T_{l}| + 1 = ({T_val}) + 1 = {U_val}]")

# Please specify a positive integer value for l to solve for.
try:
    l_value_str = input("Enter a positive integer value for l: ")
    l_value = int(l_value_str)
    if l_value > 0:
        solve_for_l(l_value)
    else:
        print("Input must be a positive integer.")
except (ValueError, TypeError):
    print("Invalid input. Please enter a positive integer.")