import itertools

def poly_to_int(p):
    """Converts a list of coefficients (p[i] is coeff of x^i) to an integer."""
    val = 0
    for i, c in enumerate(p):
        if c == 1:
            val |= (1 << i)
    return val

def int_to_poly(val):
    """Converts an integer to a list of coefficients."""
    p = []
    while val > 0:
        p.append(val & 1)
        val >>= 1
    return p if p else [0]

def poly_deg(p):
    """Returns the degree of a polynomial."""
    for i in range(len(p) - 1, -1, -1):
        if p[i] == 1:
            return i
    return -1

def poly_add(p1, p2):
    """Adds two polynomials over F_2."""
    n = max(len(p1), len(p2))
    res = [0] * n
    for i in range(n):
        c1 = p1[i] if i < len(p1) else 0
        c2 = p2[i] if i < len(p2) else 0
        res[i] = c1 ^ c2
    return res

def poly_mul(p1, p2):
    """Multiplies two polynomials over F_2."""
    if (not p1 or all(c == 0 for c in p1)) or (not p2 or all(c == 0 for c in p2)):
        return [0]
    n1 = len(p1)
    n2 = len(p2)
    res = [0] * (n1 + n2 - 1)
    for i in range(n1):
        if p1[i] == 1:
            for j in range(n2):
                if p2[j] == 1:
                    res[i+j] ^= 1
    return res

def solve():
    """Finds the least degree unit."""
    # Special polynomials
    x = [0, 1]
    x_plus_1 = [1, 1]
    x4 = [0, 0, 0, 0, 1]

    # Max degree to search for a and b
    max_deg_a = 5
    max_deg_b = 5

    for unit_deg in range(1, 10):
        # Iterate over possible degrees of a and b for this unit_deg
        max_da = unit_deg
        max_db = unit_deg - 1

        # Iterate over all polynomials a and b up to max_da and max_db
        for da in range(-1, max_da + 1):
            for db in range(-1, max_db + 1):
                # We need max(da, db+1) == unit_deg
                if max(da, db + 1) != unit_deg:
                    continue
                
                # Check for u != 1, which means b != 0 or a != 1
                if db == -1 and da <= 0 : # Excludes u=0, u=1
                    continue
                
                num_a = 2**(da + 1) if da >= 0 else 1
                num_b = 2**(db + 1) if db >= 0 else 1

                for i_a in range(num_a):
                    a_int = i_a if da >= 0 else 0
                    a = int_to_poly(a_int)
                    if poly_deg(a) != da: continue

                    for i_b in range(num_b):
                        b_int = i_b if db >= 0 else 0
                        b = int_to_poly(b_int)
                        if poly_deg(b) != db: continue
                        
                        # Norm = a^2 + a*b*x^4 + b^2*(x+1)
                        a_squared = poly_mul(a, a)
                        b_squared = poly_mul(b, b)
                        
                        term1 = a_squared
                        term2 = poly_mul(poly_mul(a, b), x4)
                        term3 = poly_mul(b_squared, x_plus_1)

                        norm = poly_add(poly_add(term1, term2), term3)
                        
                        # Norm must be 1
                        if norm == [1]:
                            a_str = "+".join([f"x^{i}" for i, c in enumerate(a) if c==1]).replace("x^0","1").replace("x^1","x")
                            if not a_str: a_str = "0"
                            b_str = "+".join([f"x^{i}" for i, c in enumerate(b) if c==1]).replace("x^0","1").replace("x^1","x")
                            if not b_str: b_str = "0"
                            print(f"Found a unit of degree {unit_deg}:")
                            print(f"u = ({a_str}) + ({b_str})y")
                            return

solve()
print("<<<3>>>")