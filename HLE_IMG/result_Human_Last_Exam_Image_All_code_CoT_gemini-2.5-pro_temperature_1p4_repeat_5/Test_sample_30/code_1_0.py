import collections

def analyze_knight_puzzles():
    """
    Solves the Knights Puzzle for the five given configurations based on a parity argument.

    The key insight is that the knight's move graph on the 4x3 board is bipartite.
    For the puzzle to be solvable (swapping k white knights with k black knights),
    two conditions must be met:
    1. The total number of moves, M, must be even. This is because the set of occupied
       squares is the same in the start and end states, and each move flips the
       parity of the number of knights on one of the partitions of the bipartite graph.
    2. The permutation of swapping k pairs of knights has a sign of (-1)^k. This must
       be achievable by M moves, leading to k_pairs = M (mod 2).

    Combining these, k_pairs must be even for a solution to exist.
    """

    # We just need to count the number of knights of one color for each configuration.
    # This number is 'k' in our argument.
    # A: 4 black, 4 white -> k = 4
    # B: 3 black, 3 white -> k = 3
    # C: 2 black, 2 white -> k = 2
    # D: 3 black, 3 white -> k = 3
    # E: 3 black, 3 white -> k = 3
    configurations = {
        "A": 4,
        "B": 3,
        "C": 2,
        "D": 3,
        "E": 3,
    }

    solvable_configs = []
    print("Analyzing solvability of Knight Puzzles A-E:")
    for name, k_pairs in configurations.items():
        # The puzzle is solvable if and only if the number of knight pairs (k) is even.
        if k_pairs % 2 == 0:
            solvability = "Solvable"
            solvable_configs.append(name)
        else:
            solvability = "Unsolvable"
        
        print(f"Configuration {name}: Has {k_pairs} pairs of knights. k is {'even' if k_pairs % 2 == 0 else 'odd'}. -> {solvability}")

    print("\nConclusion:")
    if solvable_configs:
        print(f"The solvable configurations are: {', '.join(solvable_configs)}")
    else:
        print("None of the configurations are solvable.")

# Run the analysis
analyze_knight_puzzles()
