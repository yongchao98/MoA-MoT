def solve_max_points():
    """
    This function finds the maximum number of points n by checking the logical
    constraints derived from the problem's geometric conditions.
    """

    # Let R, G, Y be the number of red, green, and yellow points.
    # The problem conditions lead to a set of logical rules:
    # 1. If R >= 3, then G must be >= 3.
    # 2. If G >= 3, then Y must be >= 3.
    # 3. If Y >= 3, then R must be >= 3.
    # 4. The case where R>=3, G>=3, and Y>=3 is impossible due to a
    #    convex hull inclusion paradox (Area(C_R) < Area(C_G) < Area(C_Y) < Area(C_R)).
    #
    # These rules together imply that we must have R<=2, G<=2, and Y<=2.
    # This program will find the maximum n by searching for the highest
    # combination of (R, G, Y) that satisfies these rules.

    max_n = 0
    best_combo = (0, 0, 0)
    search_limit = 10  # A reasonable upper bound for the search

    for R in range(search_limit):
        for G in range(search_limit):
            for Y in range(search_limit):
                # Rule out the impossible case where all are >= 3
                if R >= 3 and G >= 3 and Y >= 3:
                    continue

                # Rule out cases that violate the chain of implications
                if R >= 3 and G < 3:
                    continue
                if G >= 3 and Y < 3:
                    continue
                if Y >= 3 and R < 3:
                    continue

                # If the combination is possible, check if it's a new maximum
                current_n = R + G + Y
                if current_n > max_n:
                    max_n = current_n
                    best_combo = (R, G, Y)

    print("The analysis shows that the number of points of each color (R, G, Y) cannot be simultaneously 3 or more.")
    print("This leads to the conclusion that each color can have at most 2 points.")
    print("\nThe maximum value of n is found with the combination:")
    r_max, g_max, y_max = best_combo
    print(f"{r_max} + {g_max} + {y_max} = {max_n}")

solve_max_points()