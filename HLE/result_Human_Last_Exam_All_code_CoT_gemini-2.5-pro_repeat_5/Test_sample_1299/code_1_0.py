import math

def calculate_cardinalities(l):
    """
    Calculates the cardinalities of U_l and T_l for a given positive integer l.
    
    This function demonstrates the solution based on the derived formulas.
    The final output includes the calculated values and the expressions with numbers
    substituted in, as requested by the prompt.
    """
    
    if not isinstance(l, int) or l <= 0:
        print("Error: l must be a positive integer.")
        return

    # --- Prime Factorization ---
    # Find the exponents e_i in the prime factorization l = p1^e1 * p2^e2 * ...
    exponents = []
    n = l
    p = 2
    while p * p <= n:
        if n % p == 0:
            count = 0
            while n % p == 0:
                count += 1
                n //= p
            exponents.append(count)
        p += 1
    if n > 1:
        exponents.append(1)

    # --- Part A: |U_l| ---
    # The cardinality |U_l| is tau(l), the number of divisors of l.
    # Formula: |U_l| = (e1+1) * (e2+1) * ...
    val_A = 1
    expr_A_parts = []
    
    if l == 1:
        # For l=1, exponents list is empty. tau(1) = 1.
        expr_A = "1"
        val_A = 1
    else:
        for e in exponents:
            val_A *= (e + 1)
            expr_A_parts.append(f"({e}+1)")
        expr_A = " * ".join(expr_A_parts)

    print("A) Formula for |U_l| = tau(l) = (e_1+1)(e_2+1)...(e_s+1)")
    print(f"   For l = {l}, the calculation is:")
    print(f"   |U_{l}| = {expr_A} = {val_A}")

    print("\n----------------------------------------\n")

    # --- Part B: |T_l| ---
    # The cardinality |T_l| is calculated from its definition.
    # The derived formula is piecewise:
    # - For l = 1, |T_1| = 1
    # - For l > 1, |T_l| = ( (1+2*e1)*(1+2*e2)*... ) - 1
    # Note: This correct formula cannot be expressed simply using the suggested 'd' variable.

    if l == 1:
        val_B = 1
        expr_B = "1"
    else:
        prod_term = 1
        expr_B_parts = []
        for e in exponents:
            prod_term *= (1 + 2 * e)
            expr_B_parts.append(f"(1+2*{e})")
        expr_B = f"({' * '.join(expr_B_parts)}) - 1"
        val_B = prod_term - 1
        
    print("B) Formula for |T_l|:")
    print(f"   For l = {l}, the calculation is:")
    print(f"   |T_{l}| = {expr_B} = {val_B}")


# Since the problem asks for a general solution, we demonstrate it
# with a sample value for l, for instance l = 60 = 2^2 * 3^1 * 5^1.
# The user can change this value to solve for any other l.
example_l = 60
calculate_cardinalities(example_l)