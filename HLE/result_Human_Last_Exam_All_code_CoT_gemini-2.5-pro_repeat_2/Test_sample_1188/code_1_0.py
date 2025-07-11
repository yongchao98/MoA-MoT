def solve_ramification_jumps():
    """
    Solves the Diophantine equation for the ramification jumps.
    The equation is: 30 = 7*(t1 + 1) + 3*(t2 - t1) + (t3 - t2)
    with the constraint 1 <= t1 < t2 < t3.
    """
    print("Searching for integer solutions for ramification jumps t1, t2, t3...")
    solutions = []
    # Set a reasonable search limit for the jumps.
    limit = 30
    for t1 in range(1, limit):
        for t2 in range(t1 + 1, limit):
            for t3 in range(t2 + 1, limit):
                if 7 * (t1 + 1) + 3 * (t2 - t1) + (t3 - t2) == 30:
                    solutions.append((t1, t2, t3))

    if solutions:
        for sol in solutions:
            t1, t2, t3 = sol
            print(f"Found a solution for jumps (t1, t2, t3): ({t1}, {t2}, {t3})")
            
            # The lower filtration becomes trivial at t = t3 + 1
            t_final = t3 + 1
            print(f"The group G_i has order 8 for 0 <= i <= {t1}.")
            print(f"The group G_i has order 4 for {t1+1} <= i <= {t2}.")
            print(f"The group G_i has order 2 for {t2+1} <= i <= {t3}.")
            print(f"The group G_i is trivial for i >= {t_final}.")
            print(f"The equation is: 7 * ({t1} + 1) + 3 * ({t2} - {t1}) + ({t3} - {t2}) = {7 * (t1 + 1)} + {3 * (t2 - t1)} + {t3 - t2} = 30")
            print(f"The smallest integer t for which the filtration is trivial is t3 + 1 = {t_final}.")
            
    else:
        print("No solution found within the search limit.")

solve_ramification_jumps()
