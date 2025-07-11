def solve():
    """
    This function explains the derivation of the maximum value of n.
    """
    # Based on the derivation, at least two color sets must have size less than 3.
    # To maximize n = R + G + Y, we set the size of these two sets to 2.
    G = 2
    Y = 2

    # The problem then reduces to finding the maximum R such that any triangle
    # of red points can be covered by 2 green points.
    # This maximum value is known to be 5.
    R = 5

    # The maximum value of n is the sum of these counts.
    n_max = R + G + Y

    # Output the final equation with the numbers.
    print(f"The maximum value of n is found from the optimal distribution of points.")
    print(f"The analysis shows that at least two of the color counts must be less than 3.")
    print(f"To maximize n, we set G = {G} and Y = {Y}.")
    print(f"The maximum number of red points R that can be accommodated is {R}.")
    print(f"Final Equation: n_max = R + G + Y = {R} + {G} + {Y} = {n_max}")

solve()