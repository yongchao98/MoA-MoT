import math

def get_coeffs(p, q):
    """Computes the continued fraction coefficients for p/q."""
    coeffs = []
    if q == 0:
        return []
    if q < 0 or p < 0:
         return []
    while q > 0:
        a = p // q
        coeffs.append(a)
        p, q = q, p % q
    return coeffs

def solve():
    """
    This program counts the number of 2-bridge links K(p/q) with an even numerator 'p'
    and a crossing number of at most 13. This corresponds to the topological property
    of admitting two disjoint, non-parallel, minimal genus Seifert surfaces.
    """
    max_crossing_number = 13
    p_limit = 400  # A safe upper bound for p for crossing number <= 13
    found_links = set()

    # Iterate through even p
    for p in range(2, p_limit + 1, 2):
        # Iterate through q
        for q in range(1, p):
            if math.gcd(p, q) == 1:
                coeffs = get_coeffs(p, q)
                crossing_number = sum(coeffs)

                if crossing_number <= max_crossing_number:
                    # Found a qualifying link. Now find its canonical representation
                    # to count only unique links. A link K(p,q) and its mirror
                    # K(p,p-q) are considered the same. Also, K(p,q) is the same
                    # as K(p, q_inv).
                    try:
                        q_inv = pow(q, -1, p)
                        
                        # The set of equivalent q values for an unoriented link
                        q_equivalents = {q, p - q, q_inv, p - q_inv}
                        q_canonical = min(q_equivalents)
                        
                        found_links.add((p, q_canonical))
                    except ValueError:
                        # Should not happen since gcd(p,q)=1
                        pass

    print("The 2-bridge links K(p,q) with p-even and c<=13 are counted below.")
    print("Each '1' in the sum represents a unique link found:")
    
    sorted_links = sorted(list(found_links))
    if sorted_links:
        equation_str = " + ".join(["1"] * len(sorted_links))
        print(f"{equation_str} = {len(sorted_links)}")
    else:
        print("0")

    # Uncomment the following lines to see the list of canonical K(p,q) pairs
    # print("\nCanonical representations K(p,q) of the found links:")
    # for p_val, q_val in sorted_links:
    #     print(f"K({p_val}, {q_val})")

solve()
<<<151>>>