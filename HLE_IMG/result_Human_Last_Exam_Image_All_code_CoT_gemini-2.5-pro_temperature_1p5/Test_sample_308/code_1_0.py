def solve_puzzle():
    """
    This function determines the matching between roller configurations (1-8)
    and displacement plots (A-H) and prints the resulting sequence.
    """
    # Based on the detailed step-by-step analysis:
    # A -> 7 (Irregular 4-feature driver -> non-uniform 4-cycle plot)
    # B -> 5 (6-lobe driver -> 6-cycle plot)
    # C -> 3 (2-lobe driver, high variation -> 2-cycle high-amplitude plot)
    # D -> 4 (4-lobe driver, high variation -> 4-cycle high-amplitude plot)
    # E -> 6 (4-lobe driver, smooth -> 4-cycle smooth plot)
    # F -> 8 (3-lobe driver -> 3-cycle plot)
    # G -> 1 (2-fold asymmetric driver -> 2-cycle low-amplitude plot)
    # H -> 2 (This is the tricky part: 6-lobe driver -> 4-cycle plot)

    # The mapping from plot letter to configuration number
    mapping = {
        'A': 7,
        'B': 5,
        'C': 3,
        'D': 4,
        'E': 6,
        'F': 8,
        'G': 1,
        'H': 2
    }

    # Generate the sequence of numbers for plots A through H
    result_sequence = ""
    plot_letters = sorted(mapping.keys())
    for letter in plot_letters:
        result_sequence += str(mapping[letter])

    print("The sequence of configuration numbers for plots A-H is:")
    print(result_sequence)

solve_puzzle()