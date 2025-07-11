import math

def solve_polynomial_problem():
    """
    Solves the problem of finding max(A)^min(A) - |A| for the given polynomial.
    """
    F_ORDER = 7

    def poly_div(N, D, m):
        """
        Performs polynomial division N(x) / D(x) in F_m[x].
        N, D are lists of coefficients (highest degree first).
        Returns the remainder as a list of coefficients.
        """
        # Ensure coefficients are in the field
        N = [c % m for c in N]
        D = [c % m for c in D]
        
        degN = len(N) - 1
        degD = len(D) - 1
        
        if degN < degD:
            return N
            
        # We only use monic divisors, so D[0]=1 and its inverse is 1.
        R = list(N)
        
        for i in range(degN - degD + 1):
            coeff = R[i]
            if coeff == 0:
                continue
            
            # Subtract coeff * D from R, shifted
            for j in range(1, degD + 1):
                R[i+j] = (R[i+j] - coeff * D[j]) % m
        
        # The remainder is the last degD coefficients of R
        rem = R[-(degD):]
        # Trim leading zeros from remainder
        first_nonzero = -1
        for i in range(len(rem)):
            if rem[i] != 0:
                first_nonzero = i
                break
        
        if first_nonzero == -1:
            return [0]
        else:
            return rem[first_nonzero:]

    # Step 1: Generate all irreducible monic quadratics x^2 + bx + c over F_7
    irred_quads = []
    for b in range(F_ORDER):
        for c in range(F_ORDER):
            is_reducible = False
            # Check for roots
            for x in range(F_ORDER):
                if (x*x + b*x + c) % F_ORDER == 0:
                    is_reducible = True
                    break
            if not is_reducible:
                irred_quads.append([1, b, c])

    set_A = []
    # Step 2: Iterate through all possible 'a' in F_7
    for a in range(F_ORDER):
        poly_p = [1, 0, 0, 0, a, 3]
        is_irreducible = True

        # Step 3: Check for linear factors (roots)
        for x in range(F_ORDER):
            if (pow(x, 5, F_ORDER) + a * x + 3) % F_ORDER == 0:
                is_irreducible = False
                break
        
        if not is_irreducible:
            continue

        # Step 4: Check for irreducible quadratic factors
        for quad_factor in irred_quads:
            remainder = poly_div(poly_p, quad_factor, F_ORDER)
            if remainder == [0]:
                is_irreducible = False
                break
        
        if is_irreducible:
            set_A.append(a)

    # Step 5: Perform the final calculation
    if not set_A:
        print("Set A is empty, cannot perform calculation.")
        return

    max_A = max(set_A)
    min_A = min(set_A)
    size_A = len(set_A)
    
    result = int(math.pow(max_A, min_A)) - size_A

    # Step 6: Print the results
    print(f"The set A of 'a' values for which x^5+ax+3 is irreducible is: {set_A}")
    print(f"max(A) = {max_A}")
    print(f"min(A) = {min_A}")
    print(f"|A| = {size_A}")
    print("\nThe final calculation is max(A)^min(A) - |A|:")
    print(f"{max_A}^{min_A} - {size_A} = {result}")

solve_polynomial_problem()
<<<3>>>