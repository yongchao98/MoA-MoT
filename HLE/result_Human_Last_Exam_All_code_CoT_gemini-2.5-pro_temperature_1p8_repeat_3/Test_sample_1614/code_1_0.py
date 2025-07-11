def solve_tiling_problem():
    """
    Calculates T_n, the number of ways to tile a 2xn board, for n=4.
    """
    # A dictionary to store the values of T_n.
    # T.get(k, 0) will return T[k] if k is in the dictionary, otherwise 0.
    T = {0: 1}

    # Calculate T_1
    t1 = T.get(0, 0) + T.get(-1, 0) + T.get(-3, 0)
    T[1] = t1

    # Calculate T_2
    t2 = T.get(1, 0) + T.get(0, 0) + T.get(-2, 0)
    T[2] = t2

    # Calculate T_3
    t3 = T.get(2, 0) + T.get(1, 0) + T.get(-1, 0)
    T[3] = t3

    # Calculate T_4
    t4 = T.get(3, 0) + T.get(2, 0) + T.get(0, 0)
    T[4] = t4
    
    # Print the explanation and the final calculation for T_4
    print("Let T_n be the number of ways to tile a 2xn board.")
    print("The recurrence relation is: T_n = T_{n-1} + T_{n-2} + T_{n-4}")
    print("With base cases T_0 = 1 and T_n = 0 for n < 0.")
    print("\nCalculating the values iteratively:")
    print(f"T_1 = T_0 = 1")
    print(f"T_2 = T_1 + T_0 = {T[1]} + {T[0]} = {T[2]}")
    print(f"T_3 = T_2 + T_1 = {T[2]} + {T[1]} = {T[3]}")
    print("\nTo calculate T_4:")
    print(f"T_4 = T_3 + T_2 + T_0")
    print(f"T_4 = {T[3]} + {T[2]} + {T[0]} = {T[4]}")

solve_tiling_problem()