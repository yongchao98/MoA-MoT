def solve_knights_puzzle():
    """
    Determines which of the five knight puzzle configurations are solvable.

    The solvability of this puzzle depends on a parity argument.
    A solution exists only if the number of knights of each color, k, is even.
    This is because any sequence of allowed moves corresponds to an even permutation
    of the pieces on the board, while swapping the two sets of k knights
    corresponds to a permutation with sign (-1)^k. For these to be equal,
    k must be even.
    """

    # The number of knights of each color (k) for each configuration A-E.
    # We can determine this by counting the knights in the image.
    configurations = {
        'A': 4,  # 4 white, 4 black
        'B': 4,  # 4 white, 4 black
        'C': 2,  # 2 white, 2 black
        'D': 3,  # 3 white, 3 black
        'E': 3,  # 3 white, 3 black
    }

    solvable_configs = []
    print("Analyzing the solvability of each configuration:")
    for name, k in configurations.items():
        # Check if k is even
        if k % 2 == 0:
            is_solvable = "Solvable"
            solvable_configs.append(name)
        else:
            is_solvable = "Unsolvable"
        
        print(f"Configuration {name}: Has {k} knights of each color.")
        print(f" - The number k={k} is {'even' if k % 2 == 0 else 'odd'}.")
        print(f" - Therefore, configuration {name} is {is_solvable}.\n")

    print("Conclusion: The solvable configurations are those with an even number of knights per side.")
    print(f"The solvable configurations are: {', '.join(solvable_configs)}")

solve_knights_puzzle()