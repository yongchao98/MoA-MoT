def solve_point_problem():
    """
    This function finds the maximum number of points n = R + G + Y based on the
    problem's constraints.

    The constraints derived from the problem statement are:
    1. It's not possible for R, G, and Y to all be 3 or greater simultaneously.
    2. If R >= 3, then R <= 2*G.
    3. If G >= 3, then G <= 2*Y.
    4. If Y >= 3, then Y <= 2*R.

    The code searches for the combination of (R, G, Y) that satisfies these
    rules and maximizes their sum.
    """
    max_n = 0
    best_config = (0, 0, 0)
    
    # We can limit the search space as the solution is expected to be small.
    # A limit of 30 is more than sufficient.
    search_limit = 30 
    
    for r in range(search_limit):
        for g in range(search_limit):
            for y in range(search_limit):
                
                # Rule 1: Not all can be >= 3.
                if r >= 3 and g >= 3 and y >= 3:
                    continue
                
                # Rule 2: If R >= 3, then R <= 2*G.
                cond1_ok = (r < 3) or (r <= 2 * g)
                
                # Rule 3: If G >= 3, then G <= 2*Y.
                cond2_ok = (g < 3) or (g <= 2 * y)
                
                # Rule 4: If Y >= 3, then Y <= 2*R.
                cond3_ok = (y < 3) or (y <= 2 * r)

                # If all rules are satisfied, check if we found a new maximum n.
                if cond1_ok and cond2_ok and cond3_ok:
                    current_n = r + g + y
                    if current_n > max_n:
                        max_n = current_n
                        best_config = (r, g, y)

    # Print the result as a final equation.
    r, g, y = best_config
    print(f"The maximum value of n is obtained with a configuration of (R,G,Y) = {best_config} (or its permutations).")
    print(f"The equation for the maximum n is:")
    print(f"{r} + {g} + {y} = {max_n}")

solve_point_problem()