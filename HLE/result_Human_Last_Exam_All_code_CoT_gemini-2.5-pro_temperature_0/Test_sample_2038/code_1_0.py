import math

def get_cf_sum(p, q):
    """Calculates the sum of the coefficients of the continued fraction of p/q."""
    s = 0
    while q != 0:
        s += p // q
        p, q = q, p % q
    return s

def count_knots(max_crossing):
    """
    Counts the number of 2-bridge knots with trivial Alexander polynomial
    and crossing number at most max_crossing.
    """
    knots = {}
    max_k = int(math.sqrt(max_crossing * max_crossing)) # Heuristic limit for p
    p_values = [k*k for k in range(3, max_k + 2, 2)]

    print("Finding 2-bridge knots with trivial Alexander polynomial and crossing number <= 13...")
    print("The condition for trivial Alexander polynomial is that p in K(p/q) is an odd square.\n")

    total_count = 0
    for p in p_values:
        count_for_p = 0
        knots_for_p = {}
        # Iterate through all possible q values
        for q in range(1, p):
            # A 2-bridge link K(p,q) is a knot if and only if gcd(p,q) = 1
            if math.gcd(p, q) == 1:
                c = get_cf_sum(p, q)
                if c <= max_crossing:
                    # Find the canonical representative for the knot (and its mirror)
                    try:
                        q_inv = pow(q, -1, p)
                    except ValueError:
                        # This should not happen if gcd(p,q)==1
                        continue
                    
                    # The set {q, q_inv, p-q, p-q_inv} represents the same unoriented knot
                    q_canon = min(q, q_inv, p - q, p - q_inv)
                    key = (p, q_canon)

                    # Store the knot with its minimal crossing number found so far
                    if key not in knots_for_p or c < knots_for_p[key][0]:
                        knots_for_p[key] = (c, f"K({p}/{q})")

        if knots_for_p:
            print(f"For p = {p}:")
            # Sort knots by crossing number, then by q
            sorted_knots = sorted(knots_for_p.items(), key=lambda item: (item[1][0], item[0][1]))
            for key, (c, name) in sorted_knots:
                print(f"  - Knot K({key[0]},{key[1]}) has crossing number {c} (e.g., from {name})")
            count_for_p = len(knots_for_p)
            print(f"Found {count_for_p} distinct knots for p = {p}.\n")
            total_count += count_for_p
            knots.update(knots_for_p)

    print("---")
    print(f"Total number of such knots = {len(knots)}")
    return len(knots)

# Run the calculation for crossing number at most 13
final_count = count_knots(13)
print(f"\nFinal Answer: The total number of such knots is {final_count}.")
print(f"<<<{final_count}>>>")
