def solve_caesar_cipher_puzzle():
    """
    Calculates the length of the longest possible message for Caesar's new cipher.
    """
    # The total number of characters the paper can hold. This is the maximum length
    # of the final, encrypted message.
    paper_capacity = 10000

    # To find the length of the *longest possible* original message, we must use an
    # encryption that is as compact as possible. This means each character in the
    # original message should be represented by the minimum number of Roman numerals.
    # The absolute minimum is 1. For example, 'A' -> 'I', 'B' -> 'V', etc.
    min_encrypted_chars_per_original_char = 1

    # The maximum length of the original message is the total paper capacity
    # divided by the number of encrypted characters used for each original character.
    max_message_length = paper_capacity // min_encrypted_chars_per_original_char

    # We print the final equation showing each number.
    print(f"The maximum length of the encrypted message is {paper_capacity} characters.")
    print(f"The most efficient encryption uses {min_encrypted_chars_per_original_char} Roman numeral per original letter.")
    print(f"Therefore, the calculation for the longest message is:")
    print(f"{paper_capacity} / {min_encrypted_chars_per_original_char} = {max_message_length}")
    print(f"\nThe length of his longest message is {max_message_length}.")

solve_caesar_cipher_puzzle()