def calculate_tiling_T4():
    """
    Calculates T_4 for the tiling problem using the derived recurrence relation.
    The recurrence is T_n = T_{n-1} + 2*T_{n-2} + T_{n-4}.
    """
    
    # A dictionary is used to store the values of T_n (memoization).
    T = {0: 1} # Base case T_0 = 1

    # A helper function to compute T_n using the recurrence relation.
    # It returns 0 for n < 0.
    def get_T(n):
        if n < 0:
            return 0
        if n not in T:
            T[n] = get_T(n - 1) + 2 * get_T(n - 2) + get_T(n - 4)
        return T[n]

    # Calculate T_n for n=1, 2, 3, 4.
    T_1 = get_T(1)
    T_2 = get_T(2)
    T_3 = get_T(3)
    T_4 = get_T(4)
    T_0 = get_T(0)

    # Output the final equation with the computed values.
    # This fulfills the requirement to show each number in the final equation.
    print(f"To find T_4, we use the recurrence T_n = T_(n-1) + 2*T_(n-2) + T_(n-4).")
    print(f"We need the values for T_3, T_2, and T_0.")
    print(f"T_0 = {T_0}")
    print(f"T_1 = T_0 = {T_1}")
    print(f"T_2 = T_1 + 2*T_0 = {T_1} + 2*{T_0} = {T_2}")
    print(f"T_3 = T_2 + 2*T_1 = {T_2} + 2*{T_1} = {T_3}")
    print(f"\nNow we calculate T_4:")
    print(f"T_4 = T_3 + 2 * T_2 + T_0")
    print(f"T_4 = {T_3} + 2 * {T_2} + {T_0} = {T_4}")


# Execute the function to find the solution.
calculate_tiling_T4()