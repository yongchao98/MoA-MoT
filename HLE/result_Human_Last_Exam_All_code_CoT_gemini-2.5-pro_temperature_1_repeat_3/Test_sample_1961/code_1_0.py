import sys

def solve():
    """
    This function searches for integer solutions to the equation P(p) + P(q) = P(r)
    for integers p, q, r > 1, where P(t) = 2t^3 + 3t^2 + 6t - 11.
    It seeks the solution that has the smallest possible value for r.
    """

    def P(t):
        """Calculates the polynomial value P(t) = 2t^3 + 3t^2 + 6t - 11."""
        return 2 * t**3 + 3 * t**2 + 6 * t - 11

    # We search for solutions within a reasonable integer range.
    # The problem asks for the smallest r, so we prioritize checking smaller values.
    limit = 200
    
    # Create a mapping from P(r) values to r, for efficient lookup.
    # We are interested in r > 1.
    p_values = {P(i): i for i in range(2, limit + 1)}
    
    found_solutions = []

    # Iterate through possible values for p and q (p, q > 1).
    # To avoid duplicate pairs (p,q) and (q,p), we enforce q >= p.
    for p in range(2, limit + 1):
        for q in range(p, limit + 1):
            target_pr = P(p) + P(q)
            
            # Check if this sum corresponds to a P(r) for some integer r.
            if target_pr in p_values:
                r = p_values[target_pr]
                # Ensure r > 1, which is guaranteed by how p_values is constructed.
                found_solutions.append((p, q, r))

    if not found_solutions:
        print("no")
    else:
        # Sort solutions by the value of r to find the minimum.
        found_solutions.sort(key=lambda x: x[2])
        
        # The first solution in the sorted list has the smallest r.
        p_sol, q_sol, r_sol = found_solutions[0]
        
        # Print the numbers in the final equation as requested.
        # The final equation is P(p) + P(q) = P(r).
        print(f"A solution exists. The equation is P(p) + P(q) = P(r), where P(t) = 2t^3 + 3t^2 + 6t - 11.")
        print(f"The solution with the smallest r > 1 is (p, q, r) = ({p_sol}, {q_sol}, {r_sol}).")
        print(f"P({p_sol}) + P({q_sol}) = {P(p_sol)} + {P(q_sol)} = {P(r_sol)} = P({r_sol}).")
        print(f"The smallest possible value of r is {r_sol}.")
        
        # Final answer format
        print(f"<<<{r_sol}>>>")

solve()