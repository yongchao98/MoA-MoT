import sys

def pol_add(p1, p2):
    """Adds two polynomials in F_2[x] represented as integers."""
    return p1 ^ p2

def pol_mul(p1, p2):
    """Multiplies two polynomials in F_2[x] represented as integers."""
    res = 0
    p1_copy = p1
    while p2 > 0:
        if p2 & 1:
            res ^= p1_copy
        p1_copy <<= 1
        p2 >>= 1
    return res

def pol_deg(p):
    """Calculates the degree of a polynomial represented as an integer."""
    if p == 0:
        return -1
    return p.bit_length() - 1

def pol_to_str(p):
    """Converts a polynomial represented as an integer to a string."""
    if p == 0:
        return "0"
    if p == 1:
        return "1"
    
    res = []
    d = pol_deg(p)
    for i in range(d, -1, -1):
        if (p >> i) & 1:
            if i == 0:
                res.append("1")
            elif i == 1:
                res.append("x")
            else:
                res.append(f"x^{i}")
    return " + ".join(res)

def find_least_degree_unit():
    """
    Finds the least degree of a non-trivial unit u = a(x) + b(x)y by searching
    iteratively through increasing degrees.
    """
    x = 2
    x4 = pol_mul(x, pol_mul(x, pol_mul(x, x)))
    x_plus_1 = pol_add(x, 1)

    max_search_degree = 10 

    for unit_deg in range(1, max_search_degree):
        # Iterate over possible degrees for b(x)
        # deg(b)+1 <= unit_deg  => deg(b) <= unit_deg - 1
        for db in range(-1, unit_deg):
            # Determine required degree for a(x)
            # max(da, db+1) = unit_deg
            if db + 1 == unit_deg:
                da_max = unit_deg
            else: # db + 1 < unit_deg
                da_max = unit_deg
            
            # Iterate through polynomials a(x) and b(x)
            start_a = (1 << da_max) if da_max >= 0 else 0
            end_a = (1 << (da_max + 1)) if da_max >= 0 else 1
            for a in range(start_a, end_a):
                if pol_deg(a) != da_max: continue # ensure exact degree

                start_b = (1 << db) if db >= 0 else 0
                end_b = (1 << (db + 1)) if db >= 0 else 1
                for b in range(start_b, end_b):
                    if pol_deg(b) != db: continue # ensure exact degree

                    # We are looking for u != 1
                    if a == 1 and b == 0:
                        continue
                    
                    # Check if max(deg(a), deg(b)+1) matches the target unit_deg
                    current_unit_deg = max(pol_deg(a), pol_deg(b) + 1)
                    if current_unit_deg != unit_deg:
                        continue

                    # Calculate Norm N(a+by) = a^2 + abx^4 + b^2(x+1)
                    a_sq = pol_mul(a, a)
                    b_sq = pol_mul(b, b)
                    ab = pol_mul(a, b)
                    
                    term1 = a_sq
                    term2 = pol_mul(ab, x4)
                    term3 = pol_mul(b_sq, x_plus_1)
                    
                    norm = pol_add(pol_add(term1, term2), term3)
                    
                    if norm == 1:
                        print(f"Found a unit u = a(x) + b(x)y with least degree {unit_deg}.")
                        print(f"a(x) = {pol_to_str(a)}")
                        print(f"b(x) = {pol_to_str(b)}")
                        print(f"The degree of the unit is {unit_deg}.")
                        return
    
    print("No unit found within the search limit.")

if __name__ == '__main__':
    find_least_degree_unit()
