def solve_knights_puzzle():
    """
    Determines which configurations of the Knights Puzzle are solvable.

    The solvability of this puzzle depends on a key invariant related to permutations.
    The goal is to swap k white knights with k black knights. This is a permutation of pieces
    that has a sign of (-1)^k.
    Each move is a swap between a knight and an empty square, which is a transposition. A sequence of M moves
    results in a permutation of sign (-1)^M.
    Through analysis of the board's bipartite nature and the alternating turn rule, it can be shown that
    the total number of moves M must be even.
    Therefore, for a solution to exist, we must have (-1)^k = (-1)^M = 1, which implies k must be even.

    This function checks the number of knight pairs (k) for each configuration.
    """
    configurations = {
        'A': {
            'white': [(0, 2), (1, 2), (2, 2), (3, 2)],
            'black': [(0, 0), (1, 0), (2, 0), (3, 0)],
        },
        'B': {
            'white': [(1, 1), (3, 0), (3, 2)],
            'black': [(0, 1), (2, 0), (2, 2)],
        },
        'C': {
            'white': [(0, 0), (2, 1)],
            'black': [(0, 2), (1, 2)],
        },
        'D': {
            'white': [(0, 1), (2, 1)],
            'black': [(1, 1), (3, 1)],
        },
        'E': {
            'white': [(0, 1), (0, 2), (1, 2)],
            'black': [(0, 0), (1, 0), (1, 1)],
        }
    }

    solvable_configs = []
    print("Analyzing solvability of each configuration:")
    for name, config in configurations.items():
        num_white_knights = len(config['white'])
        num_black_knights = len(config['black'])

        # The problem states the number of knights is the same.
        # We can assert this for correctness.
        assert num_white_knights == num_black_knights
        k = num_white_knights

        # The condition for solvability is that k (the number of pairs) must be even.
        is_solvable = (k % 2 == 0)

        print(f"Configuration {name}:")
        print(f"  - Number of white knights: {num_white_knights}")
        print(f"  - Number of black knights: {num_black_knights}")
        print(f"  - Number of knight pairs (k): {k}")
        if is_solvable:
            print(f"  - k ({k}) is even. Conclusion: Solvable.")
            solvable_configs.append(name)
        else:
            print(f"  - k ({k}) is odd. Conclusion: Unsolvable.")
        print("-" * 20)

    print("\nSummary:")
    print(f"The solvable configurations are: {', '.join(solvable_configs)}")


solve_knights_puzzle()