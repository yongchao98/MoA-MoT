def solve_longest_message_length():
    """
    Calculates the length of the longest message Caesar can write under the given constraints.
    """
    paper_capacity = 10000

    # To maximize the original message's length, we must use characters
    # that have the shortest encoded representation (the minimum "cost").
    # A space character has a length of 1.
    # For letters (A-Z -> 1-26), the shortest Roman numeral representations are
    # 'I' (for A), 'V' (for E), and 'X' (for J), which are all 1 character long.
    min_encoded_char_length = 1

    # The maximum message length is the total capacity divided by the minimum cost per character.
    max_message_length = paper_capacity // min_encoded_char_length

    print(f"The total character capacity of the paper is {paper_capacity}.")
    print(f"The minimum number of characters required to write a single letter or space is {min_encoded_char_length}.")
    print("By using only characters with this minimum cost, we can maximize the message length.")
    print("\nThe final equation to calculate the longest message length is:")
    print(f"{paper_capacity} / {min_encoded_char_length} = {max_message_length}")

solve_longest_message_length()
<<<10000>>>