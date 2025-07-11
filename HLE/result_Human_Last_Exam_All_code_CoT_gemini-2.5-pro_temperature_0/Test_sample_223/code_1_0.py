def solve_chair_riddle():
    """
    Solves the chair riddle by simulating the optimal strategy.
    """
    num_chairs = 20
    chairs = ['E'] * num_chairs  # 'E' for Empty

    # The optimal strategy is to fill every other chair.
    # This ensures that when each person sits, they have no neighbors.
    # This is the only way to increase the number of occupied chairs.
    for i in range(0, num_chairs, 2):
        chairs[i] = 'O'  # 'O' for Occupied

    occupied_count = chairs.count('O')

    print(f"There are {num_chairs} chairs in a row.")
    print("The optimal strategy is to place a person on every other chair.")
    print("This results in the following final configuration:")
    print(f"Final state: {chairs}")
    print("\nWith this configuration, no more people can be added to increase the total count.")
    
    # The "equation" shows the logic: 20 chairs form 10 groups of 2 (O E).
    print("\nThe calculation is based on grouping chairs:")
    group_size = 2
    num_groups = num_chairs / group_size
    print(f"Total chairs ({num_chairs}) / Chairs per group ({group_size}) = Maximum occupied chairs ({int(num_groups)})")
    
    print(f"\nTherefore, the maximum number of chairs that can be occupied is {occupied_count}.")

solve_chair_riddle()