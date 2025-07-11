import itertools

def poly_add(p1, p2):
    """Adds two polynomials in F_2."""
    res = [0] * max(len(p1), len(p2))
    for i in range(len(p1)):
        res[i] ^= p1[i]
    for i in range(len(p2)):
        res[i] ^= p2[i]
    # Trim leading zeros
    while len(res) > 1 and res[-1] == 0:
        res.pop()
    return res

def poly_mul(p1, p2):
    """Multiplies two polynomials in F_2."""
    if (len(p1) == 1 and p1[0] == 0) or (len(p2) == 1 and p2[0] == 0):
        return [0]
    res = [0] * (len(p1) + len(p2) - 1)
    for i in range(len(p1)):
        for j in range(len(p2)):
            if p1[i] == 1 and p2[j] == 1:
                res[i+j] ^= 1
    return res

def find_unit():
    """Finds the unit of least degree."""
    # Polynomials represented as lists of coefficients (lowest degree first)
    # e.g., x^2+1 is [1, 0, 1]
    X = [0, 1]
    X4 = [0, 0, 0, 0, 1]
    X_plus_1 = [1, 1]
    ONE = [1]
    
    # Iterate through possible degrees D of the unit u
    for D in range(10): # Search up to degree 9
        # Case 1: da = D, db = D+3
        da = D
        db = D + 3
        
        # Iterate over all polynomials a(x) of degree da
        # Leading coefficient of a(x) is 1
        a_coeffs = [1] * (da)
        for a_c in itertools.product([0, 1], repeat=len(a_coeffs)):
            a = list(a_c) + [1]
            
            # Iterate over all polynomials b(x) of degree db
            # Leading coefficient of b(x) is 1
            b_coeffs = [1] * (db)
            for b_c in itertools.product([0, 1], repeat=len(b_coeffs)):
                b = list(b_c) + [1]

                # Check the norm equation: a^2 + a*b*x^4 + b^2*(x+1) = 1
                a_sq = poly_mul(a, a)
                b_sq = poly_mul(b, b)
                
                term1 = a_sq
                term2 = poly_mul(poly_mul(a, b), X4)
                term3 = poly_mul(b_sq, X_plus_1)
                
                norm_val = poly_add(poly_add(term1, term2), term3)
                
                if norm_val == ONE:
                    print(f"Found a unit of degree D = {D}")
                    # The question only asks for the degree, not the unit itself.
                    return D

        # Case 2 for completeness (will find higher degree solutions first)
        if D >= 4:
            da = D
            db = D - 4
            a_coeffs = [1] * (da)
            for a_c in itertools.product([0, 1], repeat=len(a_coeffs)):
                # Handle db=0 case
                if db < 0: continue
                a = list(a_c) + [1]
                if db == 0:
                    b_list = [[1]]
                else:
                    # Iterate over b(x)
                    b_coeffs = [1] * (db)
                    b_list = (list(b_c) + [1] for b_c in itertools.product([0, 1], repeat=len(b_coeffs)))
                
                for b in b_list:
                    a_sq = poly_mul(a, a)
                    b_sq = poly_mul(b, b)
                    term1 = a_sq
                    term2 = poly_mul(poly_mul(a, b), X4)
                    term3 = poly_mul(b_sq, X_plus_1)
                    norm_val = poly_add(poly_add(term1, term2), term3)
                    if norm_val == ONE:
                        print(f"Found a unit of degree D = {D}")
                        return D
    return -1

degree = find_unit()
print(f"The least degree of a unit u != 1 is {degree}.")
