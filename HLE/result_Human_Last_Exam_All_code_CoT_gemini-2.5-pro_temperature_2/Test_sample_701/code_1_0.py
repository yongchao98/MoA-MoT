def poly_to_int(p):
    """Converts a list of coefficients into an integer."""
    res = 0
    for i, c in enumerate(p):
        if c == 1:
            res |= (1 << i)
    return res

def int_to_poly(n):
    """Converts an integer into a list of coefficients."""
    if n == 0:
        return [0]
    p = []
    while n > 0:
        p.append(n & 1)
        n >>= 1
    return p

def poly_deg(p_int):
    """Calculates the degree of a polynomial given as an integer."""
    if p_int == 0:
        return -1
    return p_int.bit_length() - 1

def poly_mul(p, q):
    """Multiplies two polynomials (integer representation) in F_2[x]."""
    res = 0
    temp_q = q
    i = 0
    while (temp_q > 0):
        if (temp_q & 1):
            res ^= (p << i)
        temp_q >>= 1
        i += 1
    return res

def find_unit():
    """Finds the unit with the least degree."""
    # Max degree of u to check
    for u_deg in range(1, 10):
        # deg(u) = max(deg(p), deg(q)+1)
        # Iterate over deg(q) from -1 up to u_deg-1
        for q_deg in range(-1, u_deg):
            p_deg = u_deg if q_deg < u_deg - 1 else u_deg-1
            
            # Iterate through all polynomials q of degree q_deg
            q_start = (1 << q_deg) if q_deg >= 0 else 0
            q_end = (1 << (q_deg + 1)) if q_deg >= 0 else 1
            
            for q_int in range(q_start, q_end):
                # Iterate through all polynomials p of degree p_deg
                p_start = (1 << p_deg) if p_deg >= 0 else 0
                p_end = (1 << (p_deg + 1)) if p_deg >= 0 else 1

                for p_int in range(p_start, p_end):
                    if p_int == 1 and q_int == 0:
                        continue # Skip trivial unit

                    p2 = poly_mul(p_int, p_int)
                    
                    x4 = (1 << 4)
                    pqx4 = poly_mul(poly_mul(p_int, q_int), x4)

                    x_plus_1 = (1 << 1) | 1
                    q2 = poly_mul(q_int, q_int)
                    q2_x_plus_1 = poly_mul(q2, x_plus_1)
                    
                    norm = p2 ^ pqx4 ^ q2_x_plus_1
                    
                    if norm == 1:
                        # Found a unit
                        unit_degree = max(poly_deg(p_int), poly_deg(q_int)+1)
                        print(f"Found a non-trivial unit of degree {unit_degree}:")
                        
                        # Pretty print the polynomials
                        p_poly_str = "0" if p_int==0 else "+".join([f"x^{i}" if i > 1 else ("x" if i==1 else "1") for i, c in enumerate(int_to_poly(p_int)) if c == 1][::-1])
                        q_poly_str = "0" if q_int==0 else "+".join([f"x^{i}" if i > 1 else ("x" if i==1 else "1") for i, c in enumerate(int_to_poly(q_int)) if c == 1][::-1])

                        print(f"u = ({p_poly_str}) + ({q_poly_str})y")
                        return unit_degree

    return -1

min_degree = find_unit()
print(f"The least degree of a unit u != 1 is: {min_degree}")

<<<3>>>