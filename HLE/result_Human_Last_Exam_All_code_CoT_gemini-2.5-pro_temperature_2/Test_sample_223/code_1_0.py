def solve_chair_puzzle():
    """
    Solves the chair puzzle by simulating the optimal strategy and explaining the result.
    """
    n_chairs = 20
    # A list to represent the chairs. 0 is empty, 1 is occupied.
    chairs = [0] * n_chairs
    occupied_count = 0

    print("Solving the puzzle for 20 chairs.")
    print("The goal is to maximize the number of occupied chairs.")
    print("The number of occupants increases only if a new person sits on a chair with no occupied neighbors.\n")
    print("Let's simulate the optimal strategy: placing people with one empty chair in between.\n")

    # The optimal strategy is to occupy chairs 1, 3, 5, ... (or 0, 2, 4, ... in 0-based index)
    for i in range(10):
        # Calculate the 0-based index of the chair to occupy
        # For step 1 (i=0), pos=0. For step 2 (i=1), pos=2, etc.
        pos_to_sit = i * 2
        
        # Before this move, chair at pos_to_sit and its neighbors were empty in our strategy
        chairs[pos_to_sit] = 1
        occupied_count += 1
        
        # We present the step number from 1 to 10
        step_number = i + 1
        # We present the chair number from 1 to 20
        chair_number = pos_to_sit + 1
        
        print(f"Step {step_number}: A person sits on chair {chair_number}.")
        print(f"Number of occupied chairs = {occupied_count}")
        print("Chairs configuration (1=occupied, 0=empty):")
        print(f"{chairs}\n")
    
    print("--- Simulation Complete ---")
    print(f"After 10 moves, we have reached {occupied_count} occupied chairs.")
    print("The occupied chairs are at positions: 1, 3, 5, 7, 9, 11, 13, 15, 17, 19.")
    print("The empty chairs are at positions: 2, 4, 6, 8, 10, 12, 14, 16, 18, 20.\n")
    
    print("Can we add more people?")
    print("Let's check if any more 'increasing moves' are possible.")
    print("An increasing move requires an empty chair with no occupied neighbors.")
    print("However, every empty chair (e.g., chair 2, 4, etc.) now has at least one occupied neighbor.")
    print("For example, empty chair 4 has neighbors 3 and 5, which are both occupied.\n")

    print("If another person sits on any empty chair:")
    print("- One of their new neighbors (who is currently sitting) will have to leave.")
    - The number of occupants will increase by 1, then decrease by 1.")
    - The total number of occupants will remain 10.\n")
    
    print("Therefore, it's not possible to have more than 10 people at the same time.")
    print("\nThe final equation shows the number of people increasing step by step:")
    equation = " + ".join(map(str, range(1, 11)))
    print(f"1 person sits... number becomes 1.")
    print(f"Another person sits... number becomes 2.")
    print("...")
    print(f"The 10th person sits... number becomes 10.")


solve_chair_puzzle()
print("\nMaximum number of chairs that can be occupied: 10")