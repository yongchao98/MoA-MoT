def solve_knights_puzzle():
    """
    Analyzes five configurations of the Knights Puzzle on a 4x3 board
    to determine which are solvable based on a coloring invariant.
    """
    
    # 1. Model the board and its coloring
    # Board squares are numbered 0-11
    # 0  1  2
    # 3  4  5
    # 6  7  8
    # 9 10 11
    # A square's color is (row + col) % 2. 0 for black, 1 for white.
    white_squares = set()
    for r in range(4):
        for c in range(3):
            if (r + c) % 2 == 1:
                white_squares.add(r * 3 + c)

    print("Board representation:")
    print("White-colored squares are:", sorted(list(white_squares)))
    print("-" * 30)

    # 2. Define the five initial configurations from the image
    # B_pos: Black knight positions, W_pos: White knight positions
    configurations = {
        'A': {
            'B_pos': {0, 3, 6, 9},
            'W_pos': {2, 5, 8, 11}
        },
        'B': {
            'B_pos': {1, 3, 6, 8},
            'W_pos': {2, 5, 9, 11}
        },
        'C': {
            'B_pos': {2, 5},
            'W_pos': {1, 7}
        },
        'D': {
            'B_pos': {4, 10},
            'W_pos': {1, 7}
        },
        'E': {
            'B_pos': {0, 3, 4},
            'W_pos': {1, 2, 5}
        }
    }

    solvable_configs = []

    # 3. Check each configuration against the invariant
    for name, config in configurations.items():
        b_pos = config['B_pos']
        w_pos = config['W_pos']

        # Count knights on white-colored squares
        black_on_white = len(b_pos.intersection(white_squares))
        white_on_white = len(w_pos.intersection(white_squares))

        # Check for solvability
        is_solvable = (black_on_white == white_on_white)
        
        print(f"Configuration {name}:")
        print(f"  Black knights on white squares: {black_on_white}")
        print(f"  White knights on white squares: {white_on_white}")
        
        if is_solvable:
            print(f"  Result: Solvable, because {black_on_white} == {white_on_white}.")
            solvable_configs.append(name)
        else:
            print(f"  Result: Unsolvable, because {black_on_white} != {white_on_white}.")
        print("-" * 30)
    
    print("\nConclusion:")
    if solvable_configs:
        print("The solvable configuration(s) are: " + ", ".join(solvable_configs))
    else:
        print("None of the configurations are solvable.")

solve_knights_puzzle()