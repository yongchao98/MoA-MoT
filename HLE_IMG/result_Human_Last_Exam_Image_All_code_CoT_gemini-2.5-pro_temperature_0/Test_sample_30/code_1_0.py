def solve_knights_puzzle():
    """
    Solves the Knights Puzzle by analyzing the board's graph structure.

    The 4x3 board for knight moves splits into two disconnected components.
    A puzzle is solvable only if the number of white and black knights
    are equal within each component. This is a necessary and sufficient condition.
    """

    # The 12 squares of the 4x3 board can be partitioned into two
    # disconnected components for knight moves. A knight can never move
    # from a square in one component to a square in the other.
    # Component 1:
    component1 = {
        (0, 0), (0, 2), (1, 0), (1, 2), (2, 1), (3, 1)
    }

    # Define the initial configurations from the image.
    # (row, col) coordinates, with (0,0) at the top-left.
    configs = {
        'A': {
            'black': [(0, 0), (1, 0), (2, 0), (3, 0)],
            'white': [(0, 2), (1, 2), (2, 2), (3, 2)]
        },
        'B': {
            'black': [(0, 1), (1, 0), (2, 2)],
            'white': [(1, 1), (3, 0), (3, 2)]
        },
        'C': {
            'black': [(0, 2), (1, 2)],
            'white': [(0, 0), (2, 1)]
        },
        'D': {
            'black': [(1, 1), (3, 1)],
            'white': [(0, 0), (2, 1)]
        },
        'E': {
            'black': [(0, 0), (1, 0), (1, 1)],
            'white': [(0, 1), (0, 2), (1, 2)]
        }
    }

    solvable_configs = []

    print("Analyzing the five configurations of the Knights Puzzle.")
    print("A configuration is solvable if and only if the number of white and black knights in each of the two disconnected components of the board is equal.")
    print("-" * 50)

    for name, config in configs.items():
        white_pos = config['white']
        black_pos = config['black']

        # Count knights in component 1
        white_in_c1 = sum(1 for pos in white_pos if pos in component1)
        black_in_c1 = sum(1 for pos in black_pos if pos in component1)

        # Since the total number of white and black knights are equal for the whole board,
        # we only need to check the balance in one component. If it's balanced in C1,
        # it must also be balanced in C2.

        print(f"Configuration {name}:")
        print(f"  Checking the condition: (Num. white knights in C1) == (Num. black knights in C1)")
        
        if white_in_c1 == black_in_c1:
            print(f"  Equation: {white_in_c1} == {black_in_c1}")
            print("  Result: The condition is met. Configuration is SOLVABLE.")
            solvable_configs.append(name)
        else:
            print(f"  Equation: {white_in_c1} == {black_in_c1}")
            print("  Result: The condition is NOT met. Configuration is UNSOLVABLE.")
        print("-" * 20)

    print("\nFinal Answer:")
    print("The solvable configurations are:", ", ".join(sorted(solvable_configs)))

# Execute the analysis
solve_knights_puzzle()
<<<A, C, E>>>