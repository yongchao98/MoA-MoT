import numpy as np

def poly_add(p1, p2):
    """Adds two polynomials in F_2[x]."""
    l1, l2 = len(p1), len(p2)
    res = [0] * max(l1, l2)
    for i in range(len(res)):
        c1 = p1[i] if i < l1 else 0
        c2 = p2[i] if i < l2 else 0
        res[i] = (c1 + c2) % 2
    # Trim leading zeros
    while len(res) > 1 and res[-1] == 0:
        res.pop()
    return tuple(res)

def poly_mul(p1, p2):
    """Multiplies two polynomials in F_2[x]."""
    if p1 == (0,) or p2 == (0,):
        return (0,)
    l1, l2 = len(p1), len(p2)
    res = [0] * (l1 + l2 - 1)
    for i in range(l1):
        for j in range(l2):
            res[i+j] = (res[i+j] + p1[i] * p2[j]) % 2
    # Trim leading zeros
    while len(res) > 1 and res[-1] == 0:
        res.pop()
    return tuple(res)

def find_least_degree_unit():
    """
    Finds the least degree of a non-trivial unit in the ring.
    """
    x = (0, 1)
    x_plus_1 = (1, 1)
    x_4 = (0, 0, 0, 0, 1)
    one = (1,)
    
    max_deg_check = 10 # Search up to degree 10 for the unit
    
    for unit_deg in range(1, max_deg_check + 1):
        # Case 1: deg(a) = deg(b) + 4
        # deg(u) = deg(a) = deg(b) + 4
        d_b = unit_deg - 4
        if d_b >= 0:
            d_a = unit_deg
            # Iterate through all polynomials b of degree d_b
            for i in range(2**(d_b + 1)):
                b = [int(c) for c in bin(i)[2:]]
                b.reverse()
                b = tuple(b)
                if len(b)-1 != d_b and d_b > -1 : continue # Ensure correct degree
                if d_b > -1 and b[d_b] == 0: continue

                # Iterate through all polynomials a of degree d_a
                for j in range(2**(d_a + 1)):
                    a = [int(c) for c in bin(j)[2:]]
                    a.reverse()
                    a = tuple(a)
                    if len(a)-1 != d_a: continue # Ensure correct degree
                    if a[d_a] == 0: continue
                    
                    a_sq = poly_mul(a, a)
                    b_sq = poly_mul(b, b)
                    
                    term1 = a_sq
                    term2 = poly_mul(poly_mul(a, b), x_4)
                    term3 = poly_mul(b_sq, x_plus_1)
                    
                    norm = poly_add(poly_add(term1, term2), term3)
                    
                    if norm == one:
                        return unit_deg

        # Case 2: deg(b) = deg(a) + 3
        # deg(u) = deg(b) + 1 = deg(a) + 4
        d_a = unit_deg - 4
        if d_a >= 0:
            d_b = unit_deg - 1
            # Iterate through all polynomials a of degree d_a
            for i in range(2**(d_a + 1)):
                a = [int(c) for c in bin(i)[2:]]
                a.reverse()
                a = tuple(a)
                if len(a)-1 != d_a and d_a > -1 : continue # Ensure correct degree
                if d_a > -1 and a[d_a] == 0: continue

                # Iterate through all polynomials b of degree d_b
                for j in range(2**(d_b + 1)):
                    b = [int(c) for c in bin(j)[2:]]
                    b.reverse()
                    b = tuple(b)
                    if len(b)-1 != d_b: continue # Ensure correct degree
                    if b[d_b] == 0: continue
                    
                    a_sq = poly_mul(a, a)
                    b_sq = poly_mul(b, b)
                    
                    term1 = a_sq
                    term2 = poly_mul(poly_mul(a, b), x_4)
                    term3 = poly_mul(b_sq, x_plus_1)
                    
                    norm = poly_add(poly_add(term1, term2), term3)
                    
                    if norm == one:
                        if a == (1,) and b == (0,): continue # Trivial unit
                        return unit_deg

    return None

result = find_least_degree_unit()
if result:
    print(f"{result}")
else:
    print("No unit found within the checked degree limit.")
