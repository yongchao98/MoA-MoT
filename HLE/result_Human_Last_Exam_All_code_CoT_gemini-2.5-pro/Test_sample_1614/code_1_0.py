def calculate_tiling_ways(n):
    """
    Calculates T_n, the number of ways to tile a 2 x n board.
    The function also prints the breakdown of the calculation for T_n.
    """
    if n < 0:
        print("n must be a non-negative integer.")
        return

    # T is a list to store the results T_k, where T[k] stores the value for T_k.
    # The list size is n+1 to store values from T_0 to T_n.
    T = [0] * (n + 1)

    # Base case: T_0 = 1 (one way to tile a 2x0 board: the empty tiling).
    if n >= 0:
        T[0] = 1

    # Calculate T_k for k from 1 to n using the recurrence relation.
    # T_k = T_{k-1} + 2*T_{k-2} + T_{k-4}
    for k in range(1, n + 1):
        term1 = T[k-1] if k - 1 >= 0 else 0
        term2 = 2 * T[k-2] if k - 2 >= 0 else 0
        term3 = T[k-4] if k - 4 >= 0 else 0
        T[k] = term1 + term2 + term3

    # Display the results as requested.
    print(f"The number of ways to tile a 2 x n board, T(n), can be found with the recurrence:")
    print(f"T(n) = T(n-1) + 2*T(n-2) + T(n-4)")
    print(f"\nTo find T({n}), we first compute the necessary preceding values:")
    for i in range(n):
        print(f"T({i}) = {T[i]}")

    print(f"\nNow we calculate T({n}) using the recurrence:")
    
    # Get the values for the final equation
    t_n_minus_1 = T[n-1] if n-1 >= 0 else 0
    t_n_minus_2 = T[n-2] if n-2 >= 0 else 0
    t_n_minus_4 = T[n-4] if n-4 >= 0 else 0
    
    # Print the equation with names and then with numbers
    print(f"T({n}) = T({n-1}) + 2 * T({n-2}) + T({n-4})")
    print(f"T({n}) = {t_n_minus_1} + 2 * {t_n_minus_2} + {t_n_minus_4}")
    
    result = T[n]
    print(f"T({n}) = {result}")

# Calculate T_4
calculate_tiling_ways(4)
<<<12>>>