def solve_chairs_puzzle():
    """
    Solves the chair puzzle by implementing the optimal strategy.
    """
    num_chairs = 20
    
    # We represent the row of chairs as a list.
    # 0 means the chair is empty, 1 means it's occupied.
    chairs = [0] * num_chairs

    # The optimal strategy is to fill every other chair. This ensures no
    # newly seated person has a neighbor, so nobody has to leave.
    # We will place people at indices 0, 2, 4, ..., 18.
    for i in range(0, num_chairs, 2):
        chairs[i] = 1

    # The final number of occupied chairs is the sum of the 1s in our list.
    max_occupied = sum(chairs)

    print("There are 20 chairs in a row.")
    print("To maximize the number of people, they must sit in a way that no two people are adjacent.")
    print("This ensures that when a new person sits, they have no neighbors to displace.\n")

    print("The final optimal arrangement of occupied (1) and empty (0) chairs is:")
    print(chairs)
    print("\nEach '1' in the list above represents an occupied chair.")
    
    # Create the equation string as requested
    occupied_chairs_str = " + ".join(["1"] * max_occupied)
    
    print(f"To find the maximum number, we sum the occupied chairs:")
    print(f"{occupied_chairs_str} = {max_occupied}")

# Execute the function to solve the puzzle
solve_chairs_puzzle()
