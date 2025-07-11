import sys

def solve_markov_equation():
    """
    This function searches for the smallest integer r > 1 that satisfies the derived equation
    from the problem statement.
    The derived Diophantine equation is J(p) + J(q) - 11 = J(r), where
    J(x) = 2x^3 + 3x^2 + 6x, and p, q, r are integers greater than 1.
    """
    
    # Define the function J(x)
    def J(x):
        return 2 * x**3 + 3 * x**2 + 6 * x

    # A cache to store computed values of J(x) to speed up the search
    j_cache = {}

    # Smallest possible r is 3, because p and q must be at least 2 and p, q < r.
    r = 3
    # We will search up to a reasonable limit for r.
    search_limit = 100 
    
    while r < search_limit:
        if r not in j_cache:
            j_cache[r] = J(r)
        j_r = j_cache[r]

        # Iterate through possible values for p
        # Since J is increasing, p must be less than r
        for p in range(2, r):
            if p not in j_cache:
                j_cache[p] = J(p)
            j_p = j_cache[p]

            # Iterate through possible values for q, with p <= q < r
            for q in range(p, r):
                if q not in j_cache:
                    j_cache[q] = J(q)
                j_q = j_cache[q]

                # Check if the equation is satisfied
                if j_p + j_q - 11 == j_r:
                    # A solution is found. We print the components of the final equation
                    # and the smallest value of r.
                    print(f"Yes, such integers exist. The smallest value for r is {r}.")
                    print(f"A solution triplet is (p, q, r) = ({p}, {q}, {r}).")
                    print(f"The integer equation is: J(p) + J(q) - 11 = J(r)")
                    print(f"Substituting the values:")
                    print(f"J({p}) + J({q}) - 11 = J({r})")
                    print(f"{j_p} + {j_q} - 11 = {j_r}")
                    print(f"{j_p + j_q - 11} = {j_r}")
                    print(f"Therefore, the limit equation lim E({p}^X_n) * lim E({q}^X_n) = lim E({r}^X_n) holds.")
                    print(f"The smallest possible value of r is {r}.")
                    # Finally, output the answer in the requested format
                    print(f"<<<{r}>>>")
                    return

        r += 1

    print("No solution found within the search limit.")

solve_markov_equation()
