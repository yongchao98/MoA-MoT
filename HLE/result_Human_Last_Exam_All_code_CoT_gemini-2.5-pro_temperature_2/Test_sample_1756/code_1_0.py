def solve():
    """
    This function identifies and prints the correct statement letters from the user's list.
    """
    # A list of the letters corresponding to the statements determined to be correct.
    correct_statement_letters = [
        'A', 'C', 'E', 'F', 'G', 'M', 'O', 'P', 'Q', 'R', 'X', 'Y'
    ]

    # The list is already in lexicographic order, but we sort it to be certain.
    correct_statement_letters.sort()

    # Print the sorted letters, separated by a comma and space.
    print(", ".join(correct_statement_letters))

solve()