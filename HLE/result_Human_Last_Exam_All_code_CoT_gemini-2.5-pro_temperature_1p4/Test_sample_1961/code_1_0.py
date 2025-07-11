def solve_markov_equation():
    """
    This function searches for the smallest integer r > 1 that satisfies the equation
    H(p) + H(q) - 11 = H(r) for some integers p, q > 1, where H(x) = 2x^3 + 3x^2 + 6x.
    It iterates through possible values of r starting from 2 and, for each r,
    searches for a valid pair (p, q).
    """

    def H(x):
        """Calculates the polynomial H(x) = 2x^3 + 3x^2 + 6x."""
        if x < 2:  # p, q, r must be > 1
            return 0
        return 2 * x**3 + 3 * x**2 + 6 * x

    r = 2
    # The search limit is set high to ensure a solution is found if one exists in a reasonable range.
    search_limit = 500
    while r < search_limit:
        target_sum = H(r) + 11

        # Iterate through possible values for p.
        # Since p <= q is assumed (due to symmetry), H(p) <= target_sum / 2.
        p = 2
        while 2 * H(p) <= target_sum:
            required_h_q = target_sum - H(p)
            
            # Now, find an integer q >= p such that H(q) = required_h_q.
            # We can use binary search because H(x) is strictly increasing for x > 1.
            q_low = p
            # Set a reasonable upper bound for q search.
            # Since H(q) ~ 2q^3 and target_sum ~ 2r^3, q will be close to r.
            q_high = r + p 
            
            q_solution = -1
            while q_low <= q_high:
                q_mid = (q_low + q_high) // 2
                h_q_mid = H(q_mid)
                
                if h_q_mid == required_h_q:
                    q_solution = q_mid
                    break
                elif h_q_mid < required_h_q:
                    q_low = q_mid + 1
                else:
                    q_high = q_mid - 1

            if q_solution != -1:
                # Found a solution (p, q, r). Since we are iterating r upwards,
                # this r is the smallest possible.
                print(f"A solution exists for integers p,q,r > 1.")
                print(f"The equation derived is 2p^3 + 3p^2 + 6p + 2q^3 + 3q^2 + 6q - 11 = 2r^3 + 3r^2 + 6r.")
                print(f"The smallest value for r is part of the triplet (p, q, r) = ({p}, {q_solution}, {r}).")
                print(f"The smallest possible value of r is {r}.")
                return r
            
            p += 1
        r += 1

    print("No solution found within the search limit.")
    return None

# Execute the search and print the smallest r.
smallest_r = solve_markov_equation()
if smallest_r is None:
    # This case should not be reached if a solution exists.
    # If it is reached, it implies no solution was found.
    print("no") # Fallback as per instructions if no solution is found.

print(f">>>{smallest_r}<<<")