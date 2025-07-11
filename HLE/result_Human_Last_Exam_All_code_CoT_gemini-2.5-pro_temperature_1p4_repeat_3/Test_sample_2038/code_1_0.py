import math

def get_partitions(n, I=1):
    """
    Generates all integer partitions of n.
    This is a recursive generator function.
    """
    yield (n,)
    for i in range(I, n // 2 + 1):
        for p in get_partitions(n - i, i):
            yield (i,) + p

def cf_to_rational(cf):
    """
    Converts a continued fraction [a1, a2, ...] to a rational number p/q.
    """
    p, q = 1, 0
    for a in reversed(cf):
        p, q = a * p + q, p
    return p, q

def rational_to_cf(p, q):
    """
    Converts a rational number p/q to its continued fraction expansion.
    """
    cf = []
    while q > 0:
        a = p // q
        cf.append(a)
        p, q = q, p % q
    return cf

def gcd(a, b):
    """
    Computes the greatest common divisor of a and b.
    """
    while b:
        a, b = b, a % b
    return a

def find_canonical_knot(p, q):
    """
    Finds the canonical representation (p, q_canon) for a knot C(p/q).
    The canonical form is chosen to uniquely represent the knot and its mirror image.
    """
    if q == 0:
        return p, q
    # Generate the set of equivalent q values based on standard equivalences for 2-bridge knots.
    # C(p,q) is equivalent to C(p,q') if q*q' = +/- 1 (mod p).
    # We consider a knot and its mirror image as the same.
    try:
        q_inv = pow(q, -1, p)
        equivalent_qs = {q, p - q, q_inv, p - q_inv}
        # The canonical q is the smallest in the set.
        return p, min(equivalent_qs)
    except ValueError:
        # This case happens if gcd(p,q) != 1, which we filter out earlier.
        return p, q


def solve_knot_problem():
    """
    Main function to solve the knot theory problem.
    """
    max_crossing_number = 13
    unique_knots = set()

    # Generate all 2-bridge knots with crossing number <= 13
    for c in range(3, max_crossing_number + 1):
        for part in get_partitions(c):
            p, q = cf_to_rational(part)

            # Knots have odd p; even p corresponds to links.
            if p % 2 == 0:
                continue
            
            # p and q must be coprime.
            if gcd(p, q) != 1:
                continue

            # Add the canonical representation to a set to ensure uniqueness.
            canonical_p, canonical_q = find_canonical_knot(p, q)
            unique_knots.add((canonical_p, canonical_q))

    # Identify knots whose Alexander polynomial is monic (p is a perfect square)
    monic_poly_knots = []
    non_monic_poly_knots = []

    for p, q in sorted(list(unique_knots)):
        isqrt_p = int(math.sqrt(p))
        if isqrt_p * isqrt_p == p:
            monic_poly_knots.append((p, q))
        else:
            non_monic_poly_knots.append((p, q))
            
    # Print the results of the calculation step-by-step
    total_knots = len(unique_knots)
    num_monic = len(monic_poly_knots)
    num_non_monic = len(non_monic_poly_knots)

    print("The property of having two disjoint non-parallel minimal genus Seifert surfaces is equivalent to the knot's Alexander polynomial being non-monic.")
    print("For a 2-bridge knot C(p,q), its Alexander polynomial is non-monic if and only if p is not a perfect square.")
    print("The problem thus reduces to counting the number of 2-bridge knots with crossing number at most 13 for which p is not a perfect square.")
    print("\n----- Calculation -----")
    
    print(f"Total number of unique 2-bridge knots with crossing number at most 13 is: {total_knots}")
    print(f"Number of these knots where p is a perfect square (monic Alexander polynomial): {num_monic}")
    print(f"The number of knots that satisfy the condition (non-monic Alexander polynomial) is the difference:")
    print(f"{total_knots} - {num_monic} = {num_non_monic}")

if __name__ == '__main__':
    solve_knot_problem()