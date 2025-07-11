def solve_chair_puzzle():
    """
    Calculates the maximum number of chairs that can be occupied based on the puzzle's rules.
    """
    num_chairs = 20
    
    # Initialize a list to represent the row of chairs, with 0 for empty and 1 for occupied.
    chairs = [0] * num_chairs
    
    # The optimal strategy is to fill every other chair.
    # This ensures that when a person sits, they have no neighbors,
    # thus maximizing the number of occupied chairs.
    for i in range(0, num_chairs, 2):
        chairs[i] = 1
        
    # Calculate the total number of occupied chairs.
    max_occupied = sum(chairs)
    
    print(f"There are {num_chairs} chairs in a row.")
    print("The optimal strategy is to occupy every other chair.")
    print("Final arrangement (1=occupied, 0=empty):")
    # The * operator unpacks the list elements for printing.
    print(*chairs)
    
    print("\nThe maximum number of chairs that can be occupied is the total number of 1s.")
    
    # To satisfy the "output each number in the final equation" requirement,
    # we print the simple division that gives the answer.
    print(f"Final Equation: {num_chairs} / 2 = {max_occupied}")

solve_chair_puzzle()
<<<10>>>