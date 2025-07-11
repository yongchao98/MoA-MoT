def solve_problem_index():
    """
    Calculates the index of the given boundary-value problem.
    """
    # Step 1: Determine the dimension of the system, n.
    # The state vector is x(t) = (x_1(t), ..., x_2024(t)), which means there are
    # 2024 differential equations in the system.
    n = 2024

    # Step 2: Determine the number of linearly independent boundary conditions, k.
    # The given boundary conditions are:
    # 1. On x_1:  x_1(T) - x_1(0) = 0
    # 2. On x_2:  x_2(T) - x_2(0) = 0
    # 3. On x_2:  5 * (x_2(T) - x_2(0)) = 0  (redundant)
    # 4. On x_2:  100 * (x_2(T) - x_2(0)) = 0 (redundant)
    # 5. On x_2: -1000 * (x_2(T) - x_2(0)) = 0 (redundant)
    # 6. On x_2024: 100 * (x_2024(T) - x_2024(0)) = 0
    #
    # Counting the unique, independent conditions, we have one for x_1, one for x_2,
    # and one for x_2024.
    k = 3

    # Step 3: Calculate the index using the formula: index = k - n.
    index = k - n

    print("Calculating the index of the boundary-value problem:")
    print(f"Dimension of the system (n): {n}")
    print(f"Number of linearly independent boundary conditions (k): {k}")
    print("\nThe index is calculated using the formula: index = k - n")
    print(f"Final Equation: index = {k} - {n}")
    print(f"Result: index = {index}")


solve_problem_index()