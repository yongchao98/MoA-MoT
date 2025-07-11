def get_cycloaddition_possibilities():
    """
    This function provides four possible [mπ+nπ] descriptions for the
    dimerization of 3-oxidopyrylium.
    """

    # List of (m, n) pairs for the [mπ+nπ] cycloadditions
    possibilities = [
        (4, 2),
        (4, 4),
        (6, 4),
        (8, 2)
    ]

    print("In the absence of a dipolarophile, 3-oxidopyrylium undergoes dimerization.")
    print("This reaction can be formally described using several [mπ+nπ] cycloaddition notations.")
    print("\nHere are four possibilities:\n")

    # Loop through the list and print each possibility
    for i, (m, n) in enumerate(possibilities, 1):
        print(f"Possibility {i}: [{m}π + {n}π]")

# Execute the function to print the possibilities
get_cycloaddition_possibilities()