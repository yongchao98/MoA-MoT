def solve_caesar_cipher():
    """
    Calculates the length of the longest possible message Caesar can write.
    """
    # The total number of characters the paper can hold.
    total_capacity = 10000

    def int_to_roman(num):
        """Converts an integer to its Roman numeral representation."""
        val_map = [
            (1000, "M"), (900, "CM"), (500, "D"), (400, "CD"),
            (100, "C"), (90, "XC"), (50, "L"), (40, "XL"),
            (10, "X"), (9, "IX"), (5, "V"), (4, "IV"), (1, "I")
        ]
        roman_numeral = []
        for val, sym in val_map:
            # Append the Roman symbol for the given value as many times as needed.
            while num >= val:
                roman_numeral.append(sym)
                num -= val
        return "".join(roman_numeral)

    # Calculate the length of the encrypted form for each letter from A (1) to Z (26).
    encrypted_lengths = []
    for i in range(1, 27):
        roman_representation = int_to_roman(i)
        encrypted_lengths.append(len(roman_representation))

    # The space character is also a valid part of the message.
    # Its most efficient representation is a single space, with a length of 1.
    space_length = 1
    encrypted_lengths.append(space_length)

    # To get the longest message, we must use characters with the shortest
    # encrypted form. We find this minimum length.
    min_cost_per_char = min(encrypted_lengths)
    
    # The longest message length is the total capacity divided by the minimum cost.
    max_message_length = total_capacity // min_cost_per_char

    print("To write the longest message, Caesar must use characters that have the shortest Roman numeral representation.")
    print(f"The minimum length of an encrypted character (e.g., 'A' -> 'I', 'E' -> 'V', or a space) is: {min_cost_per_char}")
    print("The paper has a total capacity of 10000 characters.")
    print("\nThe final equation to find the maximum message length is:")
    print(f"{total_capacity} / {min_cost_per_char} = {max_message_length}")

solve_caesar_cipher()
<<<10000>>>