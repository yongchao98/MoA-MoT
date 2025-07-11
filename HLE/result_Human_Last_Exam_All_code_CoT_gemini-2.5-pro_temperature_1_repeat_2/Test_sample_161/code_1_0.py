def solve_caesar_cipher():
    """
    Calculates the length of the longest possible message under Caesar's Roman numeral encryption.
    """

    def int_to_roman(num):
        """Converts an integer to its Roman numeral representation."""
        val_map = [
            (1000, "M"), (900, "CM"), (500, "D"), (400, "CD"),
            (100, "C"), (90, "XC"), (50, "L"), (40, "XL"),
            (10, "X"), (9, "IX"), (5, "V"), (4, "IV"), (1, "I")
        ]
        roman_num = ""
        for val, numeral in val_map:
            # Append the numeral as many times as possible
            count, num = divmod(num, val)
            roman_num += numeral * count
        return roman_num

    # The paper can hold 10000 characters of the encrypted message.
    paper_capacity = 10000

    # The message can contain capitalized letters and spaces.
    allowed_chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZ "

    min_encrypted_len = float('inf')

    # Find the character that has the shortest Roman numeral representation.
    for char in allowed_chars:
        ascii_value = ord(char)
        roman_representation = int_to_roman(ascii_value)
        current_len = len(roman_representation)
        
        if current_len < min_encrypted_len:
            min_encrypted_len = current_len

    # The longest message is the total capacity divided by the minimum length of a single character.
    max_message_length = paper_capacity // min_encrypted_len

    # Output the final equation as requested.
    print(f"{paper_capacity} / {min_encrypted_len} = {max_message_length}")

if __name__ == "__main__":
    solve_caesar_cipher()