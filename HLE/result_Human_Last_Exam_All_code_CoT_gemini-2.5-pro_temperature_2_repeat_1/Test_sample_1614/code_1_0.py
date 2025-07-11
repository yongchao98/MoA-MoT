def calculate_tiling_T4():
    """
    Calculates the number of ways to tile a 2x4 board using 2x1, 2x2, and 2x4 tiles.
    It follows a step-by-step calculation based on a recurrence relation.
    """
    
    # T_n represents the number of ways to tile a 2xn board.
    # Recurrence relation: T_n = T_{n-1} + T_{n-2} + T_{n-4}

    # Store base cases and computed values in a dictionary.
    T = {}

    # Base case T_0: There is one way to tile a 2x0 board (the empty tiling).
    T[0] = 1
    
    # Calculate T_1: A 2x1 board can only be tiled with one 2x1 tile.
    T[1] = 1
    print(f"Calculating T_1: T_1 = {T[1]}")

    # Calculate T_2: Based on tiling a 2x1 board (T_1) and a 2x0 board (T_0).
    # T_2 = T_1 + T_0
    T[2] = T[1] + T[0]
    print(f"Calculating T_2: T_2 = T_1 + T_0 = {T[1]} + {T[0]} = {T[2]}")

    # Calculate T_3: A 2x4 tile cannot be used.
    # T_3 = T_2 + T_1
    T[3] = T[2] + T[1]
    print(f"Calculating T_3: T_3 = T_2 + T_1 = {T[2]} + {T[1]} = {T[3]}")

    # Calculate T_4 using the full recurrence relation:
    # T_4 = T_3 + T_2 + T_0
    T[4] = T[3] + T[2] + T[0]
    
    print("\nCalculating the final answer T_4:")
    print(f"The recurrence relation for T_4 is: T_4 = T_3 + T_2 + T_0")
    print(f"Substituting the calculated values:")
    print(f"T_4 = {T[3]} + {T[2]} + {T[0]} = {T[4]}")
    
calculate_tiling_T4()