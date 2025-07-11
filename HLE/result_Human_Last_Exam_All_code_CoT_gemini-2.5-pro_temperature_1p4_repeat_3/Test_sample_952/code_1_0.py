def solve_diamond_problem():
    """
    Calculates the largest number of diamonds D such that any valid arrangement
    of D diamonds on a 2024x2024 grid is guaranteed to be 'movable'.
    """
    N = 2024

    # The problem asks for the largest number of diamonds D such that any valid arrangement
    # has at least one 'movable' diamond. This is equivalent to finding K_min - 1,
    # where K_min is the minimum number of diamonds required to form a 'frozen'
    # (or 'jammed') arrangement, where no diamond can be moved.

    # For an N x N grid where N is even, the minimum size of a frozen configuration
    # is known to be (N/2)^2.
    
    # Step 1: Calculate K_min.
    half_n = N // 2
    k_min = half_n * half_n
    
    # Step 2: The largest value D is K_min - 1.
    answer = k_min - 1

    # Print the explanation and the final equation.
    print(f"The grid size is N = {N}.")
    print("A 'frozen' arrangement is one where no diamond can be moved to a valid adjacent cell.")
    print("The minimum number of diamonds in a frozen arrangement on an N x N grid (with N even) is K_min = (N/2)^2.")
    print("\nCalculating the equation step-by-step:")
    print(f"K_min = ({N} / 2)^2")
    print(f"K_min = {half_n}^2")
    print(f"K_min = {k_min}")
    print("\nThe problem asks for the largest number of diamonds D for which EVERY arrangement is movable.")
    print("This value is D = K_min - 1.")
    print("\nFinal calculation:")
    print(f"D = {k_min} - 1 = {answer}")

solve_diamond_problem()
<<<1024143>>>