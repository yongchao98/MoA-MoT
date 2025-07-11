import math

def G(x):
    """
    Calculates the value of the polynomial 2x^3 + 3x^2 + 6x.
    """
    return 2 * x**3 + 3 * x**2 + 6 * x

def find_smallest_r():
    """
    Searches for the smallest integer r > 1 that satisfies the equation
    G(p) + G(q) - 11 = G(r) for some integers p, q > 1.
    """
    # We search for r starting from 2. The first one we find will be the smallest.
    r = 2
    while True:
        target = G(r) + 11
        # Search for p, q >= 2. We can assume 2 <= p <= q.
        # This implies G(p) + G(q) >= 2 * G(p), so target >= 2 * G(p).
        # We can establish a search bound for p from this.
        p_bound = int((target / 4)**(1/3.0)) + 2 # G(p) ~ 2p^3
        
        for p in range(2, p_bound + 1):
            gp = G(p)
            
            # Since p <= q, G(p) <= G(q). Thus G(p) + G(p) <= G(p) + G(q) = target.
            if 2 * gp > target:
                break
                
            gq_needed = target - gp
            
            # Now we need to find if there is an integer q >= p such that G(q) = gq_needed.
            # We can find an approximate value of q and check integers around it.
            if gq_needed <= 0:
                continue
            
            q_approx = (gq_needed / 2.0)**(1/3.0)
            
            # Check a small window of integers around the approximation.
            for q_candidate in range(max(p, int(q_approx) - 2), int(q_approx) + 3):
                if q_candidate < p:
                    continue
                if G(q_candidate) == gq_needed:
                    # Found a solution
                    q = q_candidate
                    print(f"Yes, such integers exist. A solution is found for p={p}, q={q}, r={r}.")
                    print("The equation relating the integer parameters is:")
                    print(f"G(p) + G(q) - 11 = G(r), where G(x) = 2x^3 + 3x^2 + 6x.")
                    print("Plugging in the values:")
                    print(f"G({p}) + G({q}) - 11 = {G(p)} + {G(q)} - 11 = {G(p) + G(q) - 11}")
                    print(f"G({r}) = {G(r)}")
                    print(f"The equation holds.")
                    print(f"\nThe smallest possible value of r is {r}.")
                    return r
        r += 1

find_smallest_r()