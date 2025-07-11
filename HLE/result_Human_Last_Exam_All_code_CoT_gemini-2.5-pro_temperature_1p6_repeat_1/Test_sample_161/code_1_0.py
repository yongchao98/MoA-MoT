def solve_caesar_cipher():
    """
    Calculates the length of the longest possible message for Caesar.
    """

    # The maximum number of characters the encrypted message can have.
    paper_capacity = 10000

    def int_to_roman(num):
        """Converts an integer to its Roman numeral representation."""
        # Mapping of values to Roman numerals, from largest to smallest.
        # We only need values up to 90, which is the ASCII value for 'Z'.
        val_map = [
            (90, "XC"), (50, "L"), (40, "XL"), (10, "X"),
            (9, "IX"), (5, "V"), (4, "IV"), (1, "I")
        ]
        roman_num = ""
        for val, numeral in val_map:
            # Add the numeral for each time its value fits into the number
            count = num // val
            roman_num += numeral * count
            num %= val
        return roman_num

    # The set of allowed characters for the original message.
    possible_chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZ "

    min_cost = float('inf')
    cheapest_char = ''
    
    # Find the character with the shortest Roman numeral representation.
    for char in possible_chars:
        ascii_val = ord(char)
        roman_val = int_to_roman(ascii_val)
        cost = len(roman_val)
        if cost < min_cost:
            min_cost = cost
            cheapest_char = char

    # The longest message length is the total capacity divided by the minimum cost per character.
    max_length = paper_capacity // min_cost

    print("Step 1: Determine the character with the shortest encrypted form.")
    print(f"The most 'efficient' character is '{cheapest_char}', which has an ASCII value of {ord(cheapest_char)}.")
    print(f"Its encrypted Roman numeral form is '{int_to_roman(ord(cheapest_char))}'.")
    print(f"The length of this shortest form is {min_cost} characters.")
    
    print("\nStep 2: Calculate the maximum message length.")
    print(f"The paper has a total capacity of {paper_capacity} characters.")
    print("To get the longest message, we fill it with the most efficient character.")
    
    # Output the final equation with each number.
    print("\nFinal Equation:")
    print(f"({paper_capacity} / {min_cost}) = {max_length}")
    
    print(f"\nThe length of the longest message is {max_length}.")

solve_caesar_cipher()
<<<5000>>>