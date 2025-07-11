def solve_chair_puzzle():
    """
    Solves the chair puzzle by simulating the optimal strategy.
    """
    num_chairs = 20
    # 0 represents an empty chair, 1 represents an occupied chair.
    chairs = [0] * num_chairs
    occupied_count = 0
    
    print(f"There are {num_chairs} chairs in a row, all initially empty.")
    print("The goal is to maximize the number of occupied chairs.")
    print("\nBased on the rules, the number of occupants only increases when a person")
    print("sits in a chair with no occupied neighbors. The best strategy is to fill")
    print("every other chair.\n")

    print(f"Initial state: {''.join(map(str, chairs))}")
    print(f"Initial equation: 0")
    print(f"Occupied chairs: {occupied_count}")
    print("-" * 40)
    
    step = 0
    # We place people in chairs 1, 3, 5, ... (indices 0, 2, 4, ...).
    # Each step is an "increasing move".
    for i in range(0, num_chairs, 2):
        step += 1
        previous_count = occupied_count
        
        # A person sits on an empty chair at index i.
        chairs[i] = 1
        occupied_count += 1
        
        # For this strategy, neighbors are always empty, so no one leaves.
        
        print(f"Step {step}: A person sits in chair {i + 1}.")
        # Visualizing the state where 'O' is occupied and '_' is empty.
        visual_state = "".join(['O' if c == 1 else '_' for c in chairs])
        print(f"State: {visual_state}")
        print(f"Equation: {previous_count} + 1 = {occupied_count}")
        print(f"Occupied chairs: {occupied_count}")
        print("-" * 40)
        
    print("The final configuration is reached. All empty chairs are now next to")
    print("at least one occupied chair. Any further move will not increase the count.")
    print(f"\nThe maximum number of chairs that can be occupied is {occupied_count}.")

solve_chair_puzzle()