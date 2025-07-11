import collections

def solve_knights_puzzle():
    """
    Analyzes the solvability of five Knights Puzzle configurations based on a parity argument.

    A configuration is solvable only if the number of knight pairs (k) is even.
    This is because:
    1. The knight's move graph is bipartite. This implies the total number of moves (M) must be even for the empty squares to return to their initial color distribution.
    2. The goal is to swap k white knights with k black knights, which corresponds to a permutation of sign (-1)^k.
    3. The M moves correspond to a permutation of sign (-1)^M.
    4. For a solution to exist, (-1)^k = (-1)^M.
    5. Since M must be even, (-1)^M = 1.
    6. Therefore, (-1)^k must be 1, which means k must be even.
    """
    configs = {
        'A': {
            'board': (4, 3),
            'white': {(0, 2), (1, 2), (2, 2), (3, 2)},
            'black': {(0, 0), (1, 0), (2, 0), (3, 0)},
        },
        'B': {
            'board': (4, 3),
            'white': {(1, 2), (3, 1), (3, 2)},
            'black': {(0, 1), (1, 1), (2, 0)},
        },
        'C': {
            'board': (4, 3),
            'white': {(0, 0), (2, 1)},
            'black': {(0, 2), (1, 2)},
        },
        'D': {
            'board': (4, 3),
            # Positions carefully identified from the image
            'white': {(0, 0), (2, 0), (3, 2)},
            'black': {(1, 1), (2, 1), (3, 1)},
        },
        'E': {
            'board': (3, 4), # Note the different board dimension
            'white': {(0, 2), (0, 3), (1, 2), (1, 3)},
            'black': {(0, 0), (0, 1), (1, 0), (1, 1)},
        }
    }

    solvable_configs = []
    print("Analyzing solvability of Knights Puzzle configurations...")
    print("The primary condition for solvability is that the number of knight pairs, k, must be even.")
    print("-" * 30)

    for name, config in configs.items():
        # k is the number of white knights (which equals the number of black knights)
        k = len(config['white'])
        
        print(f"Configuration {name}:")
        print(f"Number of white/black knight pairs (k) = {k}")
        
        # Check if k is even
        is_even = (k % 2 == 0)
        
        print(f"The solvability condition is k % 2 == 0.")
        print(f"Checking the condition for k = {k}: {k} % 2 = {k % 2}")
        
        if is_even:
            print(f"Result: Since k is even, configuration {name} is considered SOLVABLE.")
            solvable_configs.append(name)
        else:
            print(f"Result: Since k is odd, configuration {name} is UNSOLVABLE.")
        print("-" * 30)
        
    print("Summary:")
    print(f"The solvable configurations are: {', '.join(sorted(solvable_configs))}")

solve_knights_puzzle()