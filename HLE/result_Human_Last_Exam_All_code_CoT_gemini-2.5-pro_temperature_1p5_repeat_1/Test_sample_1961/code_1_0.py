import sys

def G(x):
    """
    Calculates the value of the function G(x) = 2x^3 + 3x^2 + 6x.
    """
    return 2 * x**3 + 3 * x**2 + 6 * x

def find_smallest_r():
    """
    Searches for integer solutions (p, q, r) > 1 to the equation
    G(p) + G(q) - 11 = G(r)
    and returns the solution with the smallest r.
    """
    # Using a cache to store computed values of G(x) for efficiency
    g_cache = {}

    def get_g(val):
        if val not in g_cache:
            g_cache[val] = G(val)
        return g_cache[val]

    # Iterate r starting from 2 to find the smallest solution
    # The upper limit for r can be increased if no solution is found.
    for r in range(2, 200):
        g_r = get_g(r)
        target = g_r + 11

        # Search for p and q.
        # Since G is increasing, p and q must be less than r.
        # We search p up to r. The loop for q starts from p (for p<=q).
        for p in range(2, r):
            g_p = get_g(p)
            
            # Since p<=q, G(p)<=G(q). So G(p)+G(p) <= G(p)+G(q) = target.
            # This provides an optimization.
            if 2 * g_p > target:
                break

            g_q_target = target - g_p
            
            # Search for q such that G(q) = g_q_target
            # We can start q search from p
            for q in range(p, r):
                g_q = get_g(q)
                if g_q == g_q_target:
                    print(f"Found a solution: p = {p}, q = {q}, r = {r}.")
                    print(f"The equation is: G({p}) + G({q}) - 11 = G({r})")
                    print(f"L({p}) * L({q}) = L({r})")
                    print("Remember that the final equation for the expectation is:")
                    print(f"lim(n->inf) E({p}^X_n) * E({q}^X_n) = lim(n->inf) E({r}^X_n)")
                    return r
                elif g_q > g_q_target:
                    # Since G is strictly increasing, no further q for this p will work.
                    break
    
    return "No solution found in the checked range."

smallest_r = find_smallest_r()
print(f"The smallest possible value of r is: {smallest_r}")

print(f"\n<<< {smallest_r} >>>")
