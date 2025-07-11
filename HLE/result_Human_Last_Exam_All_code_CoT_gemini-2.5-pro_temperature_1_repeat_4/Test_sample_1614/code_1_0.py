def calculate_t4():
    """
    Calculates the number of ways to tile a 2x4 board (T_4)
    using the derived recurrence relation T_n = T_{n-1} + T_{n-2} + T_{n-4}.
    """
    # Base cases determined by enumeration
    # T_0: 1 way (empty tiling)
    # T_1: 1 way (one 2x1 tile)
    # T_2: 2 ways (two 2x1 tiles or one 2x2 tile)
    # T_3: 3 ways (three 2x1 tiles, or one 2x1 and one 2x2 in two orders)
    T = {
        0: 1,
        1: 1,
        2: 2,
        3: 3
    }

    # Get the values for T_3, T_2, and T_0
    t3 = T[3]
    t2 = T[2]
    t0 = T[0]

    # Calculate T_4 using the recurrence relation
    t4 = t3 + t2 + t0

    # Print the step-by-step calculation
    print("To find T_4, we use the recurrence relation: T_4 = T_3 + T_2 + T_0")
    print("\nFirst, we establish the required base cases:")
    print(f"T_0 = {t0}")
    print(f"T_2 = {t2}")
    print(f"T_3 = {t3}")
    
    print("\nNow, we substitute these values into the equation:")
    print(f"T_4 = {t3} + {t2} + {t0}")
    
    print("\nFinally, we compute the result:")
    print(f"T_4 = {t4}")

calculate_t4()