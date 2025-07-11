def H(x):
    """
    Calculates the value of the helper function H(x) = 2x^3 + 3x^2 + 6x.
    """
    return 2 * x**3 + 3 * x**2 + 6 * x

def find_smallest_r():
    """
    Searches for the integer solution (p, q, r) with p, q, r > 1
    to the equation H(p) + H(q) - 11 = H(r) that minimizes r.
    """
    # A search limit for r. If no solution is found, this may need to be increased.
    limit = 100
    
    # Pre-compute H(x) values for efficiency
    h_values = {i: H(i) for i in range(1, limit + 1)}
    h_set = set(h_values.values())
    
    # Iterate through r starting from 2 to find the smallest r first.
    for r in range(2, limit + 1):
        # The target sum H(p) + H(q)
        target = h_values[r] + 11
        
        # Search for p and q such that H(p) + H(q) = target.
        # We can assume p <= q due to symmetry.
        # A reasonable optimization is that p and q won't be much larger than r.
        for p in range(2, r + 2):
            h_p = h_values[p]
            
            # Required value for H(q)
            h_q = target - h_p
            
            # Check if a valid q exists
            if h_q in h_set:
                # To find q, we can invert the h_values dictionary
                h_inv = {v: k for k, v in h_values.items()}
                q = h_inv[h_q]
                
                # Ensure q > 1 and p <= q
                if q > 1 and p <= q:
                    print(f"Found a solution that gives the smallest r: p={p}, q={q}, r={r}")
                    print(f"The equation is H({p}) + H({q}) - 11 = H({r})")
                    print(f"{H(p)} + {H(q)} - 11 = {H(r)}")
                    print(f"{H(p) + H(q) - 11} = {H(r)}")
                    print(f"The smallest possible value of r is {r}.")
                    return r
    
    print("No solution found within the search limit.")
    return None

if __name__ == '__main__':
    smallest_r = find_smallest_r()
