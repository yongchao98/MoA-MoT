def int_to_roman(num):
    """Converts an integer to a Roman numeral string."""
    val_map = [
        (1000, "M"), (900, "CM"), (500, "D"), (400, "CD"),
        (100, "C"), (90, "XC"), (50, "L"), (40, "XL"),
        (10, "X"), (9, "IX"), (5, "V"), (4, "IV"),
        (1, "I")
    ]
    roman_num = ''
    for val, numeral in val_map:
        while num >= val:
            roman_num += numeral
            num -= val
    return roman_num

def solve():
    """
    Calculates the longest possible message length based on Roman numeral encryption.
    """
    paper_capacity = 10000
    min_encrypted_length = float('inf')

    # Find the minimum length for encrypted letters A-Z (numbers 1-26)
    for i in range(1, 27):
        roman_representation = int_to_roman(i)
        length = len(roman_representation)
        if length < min_encrypted_length:
            min_encrypted_length = length

    # The most efficient representation for a space would also be a single character.
    # Therefore, the absolute minimum length is 1.
    # This is confirmed by 'A' -> 'I' (length 1), 'E' -> 'V' (length 1), etc.

    # Calculate the maximum message length
    # To write the longest message, Caesar should only use characters that encrypt
    # to the shortest possible Roman numeral (length 1).
    max_message_length = paper_capacity // min_encrypted_length

    print("The total character capacity of the paper is {}.".format(paper_capacity))
    print("To write the longest message, Caesar must use characters that result in the shortest Roman numeral.")
    print("For example, 'A' (1) -> 'I' (length 1), 'E' (5) -> 'V' (length 1), or 'J' (10) -> 'X' (length 1).")
    print("The minimum length for an encrypted character is {}.".format(min_encrypted_length))
    print("\nThe length of the longest message is calculated as:")
    print("{} / {} = {}".format(paper_capacity, min_encrypted_length, max_message_length))
    
    # Final answer in the required format
    print("\n<<<{}>>>".format(max_message_length))

solve()