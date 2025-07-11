def solve_sequence_puzzle():
    """
    Solves the puzzle by identifying the underlying sequence of words.

    The pattern is that the cryptic words represent the English words for numbers
    in ascending order: ONE, TWO, THREE, FOUR. The task is to find the decoded
    form of the fourth word in the list.
    """
    # Define the sequence of decoded words
    decoded_sequence = ["ONE", "TWO", "THREE", "FOUR", "FIVE", "SIX"]

    # The four cryptic words correspond to the first four numbers.
    # We need to decode the fourth word.
    target_index = 3  # Corresponds to the fourth item (0-indexed)

    # Get the decoded word
    answer = decoded_sequence[target_index]

    print(f"The cryptic words follow a simple numerical sequence: ONE, TWO, THREE, ...")
    print(f"The fourth word in the sequence is '{answer}'.")

solve_sequence_puzzle()