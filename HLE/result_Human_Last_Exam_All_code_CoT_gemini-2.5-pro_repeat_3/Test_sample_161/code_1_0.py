def solve_caesar_love_letter():
    """
    Calculates the maximum length of a message that can be encrypted using
    Roman numerals on a paper with a limited capacity.
    """

    # The paper can hold a maximum of 10000 characters for the encrypted message.
    paper_capacity = 10000

    # The original message alphabet has 27 characters (A-Z and space).
    # The encryption alphabet has 7 characters (I, V, X, L, C, D, M).

    # To maximize the original message length, we must find the most efficient
    # encoding, which is the one that uses the fewest Roman numerals per
    # original character.
    # The shortest possible code for an original character is a single Roman numeral.
    # For example, we could define an encoding like 'A' -> 'I'.
    min_encoding_length_per_char = 1

    # If Caesar writes a message using only characters that have an encoding
    # length of 1, the length of the original message will be equal to the
    # length of the encrypted message.

    # Let L be the length of the original message.
    # The length of the encrypted message will also be L.
    # This encrypted message must fit on the paper.
    # Therefore, L <= paper_capacity.

    # The maximum length of the original message is the paper capacity divided by
    # the minimum possible encoding length for a character.
    max_message_length = paper_capacity // min_encoding_length_per_char

    print("The problem is to find the maximum length of an original message given that its encrypted version must fit within 10000 characters.")
    print("\nThe encryption uses 7 Roman numerals, while the original message uses 27 characters (A-Z, space).")
    print("To make the original message as long as possible, we need the most compact encryption.")
    print("\nThe most compact representation for any character is a single Roman numeral (length 1).")
    print("For example, we can map 'A' -> 'I', 'B' -> 'V', and so on for 7 characters.")
    print("\nIf Caesar's message consists only of these characters, each character in his message takes up only one character on the paper.")
    print("This means the length of the original message can be equal to the length of the encrypted message.")
    print("\nGiven the paper's capacity, the maximum length is determined by the following equation:")
    print(f"{paper_capacity} / {min_encoding_length_per_char} = {max_message_length}")

solve_caesar_love_letter()
<<<10000>>>