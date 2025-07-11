def get_roman_numeral_length(n):
    """
    Calculates the length of the Roman numeral representation for a given integer n.
    This function processes the number without building the actual Roman numeral string,
    making it efficient for finding just the length.
    """
    if not 0 < n < 4000:
        # Roman numerals are typically defined in this range.
        # Characters outside this range are not considered.
        return float('inf')

    # Pairs of (value, length of Roman symbol)
    # e.g., 900 is "CM" (length 2), 1000 is "M" (length 1)
    roman_map = [
        (1000, 1), (900, 2), (500, 1), (400, 2), (100, 1), (90, 2),
        (50, 1), (40, 2), (10, 1), (9, 2), (5, 1), (4, 2), (1, 1)
    ]

    length = 0
    for val, sym_len in roman_map:
        # Count how many times the symbol is used
        count = n // val
        length += count * sym_len
        # Subtract the value from the number
        n %= val
    return length

def solve_caesar_cipher_length():
    """
    Finds the longest possible message length under Caesar's Roman numeral encryption.
    """
    paper_capacity = 10000

    # Determine the shortest possible length for an encrypted character.
    # We assume the standard A=1, B=2, ..., Z=26 mapping.
    min_encrypted_length = float('inf')
    for char_code in range(1, 27):  # Numbers 1 through 26 for A-Z
        length = get_roman_numeral_length(char_code)
        if length < min_encrypted_length:
            min_encrypted_length = length

    # To write the longest message, Caesar must only use characters that
    # encrypt to this minimum length.
    # The longest message length = paper_capacity / min_encrypted_length
    max_message_length = paper_capacity // min_encrypted_length
    
    # Print the explanation and the final equation as requested.
    print(f"The total capacity of the paper is {paper_capacity} characters.")
    print(f"To achieve the longest message, each character must have the shortest possible encrypted length.")
    print(f"The minimum length of a Roman numeral for numbers 1-26 is {min_encrypted_length}.")
    print("This corresponds to letters like 'A' (1 -> I), 'E' (5 -> V), or 'J' (10 -> X).")
    print("\nThe final calculation is:")
    print(f"{paper_capacity} / {min_encrypted_length} = {max_message_length}")

solve_caesar_cipher_length()
