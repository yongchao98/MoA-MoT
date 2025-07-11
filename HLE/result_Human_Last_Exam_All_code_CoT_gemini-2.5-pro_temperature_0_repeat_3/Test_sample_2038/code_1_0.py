import math

def get_crossing_number(p, q):
    """
    Computes the crossing number of the 2-bridge knot K(p,q)
    by summing the coefficients of its continued fraction.
    """
    if q == 0:
        return p
    
    coeffs = []
    temp_p, temp_q = p, q
    while temp_q != 0:
        a = temp_p // temp_q
        coeffs.append(a)
        temp_p, temp_q = temp_q, temp_p % temp_q
    return sum(coeffs)

def solve():
    """
    Finds and counts the 2-bridge knots with crossing number at most 13
    that have two disjoint non-parallel minimal genus Seifert surfaces.
    """
    max_crossing_number = 13
    # The largest possible determinant p for a knot with crossing number C is F_{C+1}.
    # For C=13, the max p is F_14 = 377.
    max_p = 377
    
    counts_by_c = {i: 0 for i in range(3, max_crossing_number + 1)}
    
    # p must be odd for a 2-bridge knot K(p,q)
    for p in range(3, max_p + 1, 2):
        # We only need to check q up to p/2
        for q in range(1, p // 2 + 1):
            # Condition 1: p and q must be coprime
            if math.gcd(p, q) != 1:
                continue
            
            # Condition 2: q^2 = -1 (mod p)
            if (q * q) % p == p - 1:
                # Condition 3: Crossing number must be at most 13
                cn = get_crossing_number(p, q)
                if cn <= max_crossing_number:
                    if cn in counts_by_c:
                        counts_by_c[cn] += 1

    total_knots = sum(counts_by_c.values())
    
    print("Number of 2-bridge knots with two disjoint minimal genus Seifert surfaces, by crossing number (c):")
    
    equation_parts = []
    for c in range(3, max_crossing_number + 1):
        count = counts_by_c[c]
        if count > 0:
            print(f"c = {c}: {count} knot(s)")
            equation_parts.append(f"{count}")

    print("\nFinal count expressed as a sum:")
    print(f"{' + '.join(equation_parts)} = {total_knots}")

solve()