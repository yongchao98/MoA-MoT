def solve_tiling():
    """
    Calculates the number of ways to tile a 2xn board using 2x1, 2x2, and 2x4 tiles.
    This function specifically calculates T_4 based on a recurrence relation.
    """
    n = 4
    # T stores the calculated values of T_n. We start with the base case T_0 = 1.
    T = {0: 1}

    # We loop from 1 to n to calculate T_i for each i.
    for i in range(1, n + 1):
        # The recurrence relation is T_i = T_{i-1} + T_{i-2} + T_{i-4}.
        # T.get(k, 0) is used to safely access T[k], returning 0 if k is not found
        # (which elegantly handles negative indices based on our definition T_n = 0 for n < 0).
        T[i] = T.get(i-1, 0) + T.get(i-2, 0) + T.get(i-4, 0)

    # For clarity, let's print the values leading up to T_4
    print(f"T_1 = {T[1]}")
    print(f"T_2 = {T[2]}")
    print(f"T_3 = {T[3]}")

    # Now, we get the specific values needed for the T_4 calculation.
    t4 = T[4]
    t3 = T[3]
    t2 = T[2]
    t0 = T[0]

    # Finally, print the detailed equation for T_4 as requested.
    print("\nThe calculation for T_4 is:")
    print(f"T_4 = T_3 + T_2 + T_0")
    print(f"T_4 = {t3} + {t2} + {t0} = {t4}")

solve_tiling()