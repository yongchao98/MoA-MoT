def solve_max_points():
    """
    Finds the maximum number of points n by iterating through combinations
    of red (r), green (g), and yellow (y) points that satisfy derived
    combinatorial constraints.
    """
    max_n = 0
    max_config = (0, 0, 0)

    # We check a reasonable range for r, g, y.
    for r in range(11):
        for g in range(11):
            for y in range(11):
                # Constraint A: At least one color count must be less than 3.
                # This is because r,g,y >= 3 leads to a contradiction with convex hulls.
                is_A_valid = (r < 3) or (g < 3) or (y < 3)
                if not is_A_valid:
                    continue

                # Constraint B: If one color count is 2, the others must be <= 4.
                # This is from the "dividing line" argument.
                is_B_valid = True
                if r == 2 and (g > 4 or y > 4):
                    is_B_valid = False
                if g == 2 and (r > 4 or y > 4):
                    is_B_valid = False
                if y == 2 and (r > 4 or g > 4):
                    is_B_valid = False
                if not is_B_valid:
                    continue

                # If both constraints are met, check if this is a new maximum.
                n = r + g + y
                if n > max_n:
                    max_n = n
                    max_config = (r, g, y)

    print("Based on the derived constraints, the maximum value of n is found.")
    print(f"The optimal distribution of colors is (r,g,y) = {max_config} or its permutations.")
    final_r, final_g, final_y = max_config
    final_n = max_n
    
    # The final equation output
    print("The final calculation is:")
    print(f"{final_r} + {final_g} + {final_y} = {final_n}")

solve_max_points()