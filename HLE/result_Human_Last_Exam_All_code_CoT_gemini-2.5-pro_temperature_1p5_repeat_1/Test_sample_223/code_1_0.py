import math

def solve_chair_riddle():
    """
    This function solves the riddle about the maximum number of occupied chairs.
    It simulates the optimal strategy and prints the result.
    """
    num_chairs = 20
    chairs = [0] * num_chairs  # 0 represents an empty chair, 1 for occupied

    # The optimal strategy to maximize occupants is to place each new person
    # in a chair with no occupied neighbors. This avoids the rule where a
    # neighbor has to leave. The best way to achieve this is to fill
    # every other chair.
    # We will place a person in chair 1, 3, 5, ... (indices 0, 2, 4, ...).
    for i in range(0, num_chairs, 2):
        chairs[i] = 1

    # Once this configuration is reached, any new person sitting in an empty
    # chair will necessarily have an occupied neighbor, causing that neighbor
    # to leave. The total number of occupied chairs can no longer increase.
    final_occupancy = sum(chairs)
    
    print(f"There are {num_chairs} chairs in a row.")
    print("To find the maximum number of people, we follow a strategy where each person sits down without any neighbors.")
    print("This leads to the following final arrangement (1 = occupied, 0 = empty):")
    print(chairs)
    print("\nTo find the final count, we sum the occupied chairs:")

    # Create the equation string as requested
    occupied_indices = [str(c) for c in chairs if c == 1]
    equation_string = " + ".join(occupied_indices)
    
    print(f"{equation_string} = {final_occupancy}")
    print(f"\nThe maximum number of chairs that can be occupied is {final_occupancy}.")

solve_chair_riddle()