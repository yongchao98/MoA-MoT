import math

def solve_and_print():
    """
    This function solves the problem by finding the smallest integer r > 1
    that satisfies the derived Diophantine equation for some integers p, q > 1,
    and then prints the result.

    The derivation shows that the problem is equivalent to finding integers p,q,r > 1
    satisfying: h(p) + h(q) = h(r) + 11, where h(x) = 2x^3 + 3x^2 + 6x.
    We are looking for the solution with the minimum possible value of r.
    """

    def h(x):
        """Computes the helper function h(x) = 2x^3 + 3x^2 + 6x."""
        return 2 * x**3 + 3 * x**2 + 6 * x

    # We search for the smallest r, so we iterate r upwards starting from 2.
    # For each r, we check if there exist p,q > 1 satisfying the equation.
    # From the equation, it can be proved that any solution must have p < r and q < r.
    
    r = 2
    # Use a dictionary to cache values of h(x) for performance.
    h_vals = {1: 11}
    h_inv = {11: 1}

    while True:
        h_vals[r] = h(r)
        h_inv[h_vals[r]] = r

        target = h_vals[r] + 11

        # We search for p, q in [2, r-1]. Let's assume p <= q.
        # This implies h(p) <= target / 2.
        for p in range(2, r):
            h_p = h_vals.get(p)
            if h_p is None: # Should not happen with sequential search
                h_p = h(p)
                h_vals[p] = h_p
                h_inv[h_p] = p
            
            if h_p > target / 2:
                break
                
            h_q_target = target - h_p
            
            if h_q_target in h_inv:
                q = h_inv[h_q_target]
                if q >= p and q < r: # Check if q is a valid solution
                    # We found the solution with the smallest r.
                    print(f"Yes, such integers exist.")
                    print(f"The equation relating the limits is derived to be:")
                    print(f"2p^3 + 3p^2 + 6p + 2q^3 + 3q^2 + 6q - 11 = 2r^3 + 3r^2 + 6r")
                    print(f"\nThe smallest value for r is {r}, with a corresponding solution (p, q) = ({p}, {q}).")
                    print(f"Let's verify this solution:")
                    h_p_val = h(p)
                    h_q_val = h(q)
                    h_r_val = h(r)
                    lhs = h_p_val + h_q_val
                    rhs = h_r_val + 11
                    print(f"h({p}) + h({q}) = {h_p_val} + {h_q_val} = {lhs}")
                    print(f"h({r}) + h(1) = {h_r_val} + 11 = {rhs}")
                    
                    if lhs == rhs:
                        print("The equation holds.")
                    else:
                        # This case should not be reached if logic is correct
                        print("Error: The equation does not hold. Re-check logic.")
                    
                    # Return the value of r as the final answer
                    return r

        r += 1
        # A safety break to prevent an infinite loop in case of error.
        if r > 200:
            print("No solution found within the search limit.")
            return None

# Run the solver function
smallest_r = solve_and_print()
# The final answer needs to be in <<<>>> format.
# print(f"<<<{smallest_r}>>>")