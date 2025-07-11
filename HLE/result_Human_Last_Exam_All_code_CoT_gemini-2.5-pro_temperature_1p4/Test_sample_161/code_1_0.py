def main():
    """
    Calculates the longest possible message length based on a Roman numeral
    encryption scheme.
    """

    def int_to_roman(num):
        """Converts an integer to a Roman numeral string."""
        # Standard values and their Roman numeral representations
        val_map = [
            (1000, "M"), (900, "CM"), (500, "D"), (400, "CD"),
            (100, "C"), (90, "XC"), (50, "L"), (40, "XL"),
            (10, "X"), (9, "IX"), (5, "V"), (4, "IV"),
            (1, "I")
        ]
        roman_str = ""
        for val, numeral in val_map:
            while num >= val:
                roman_str += numeral
                num -= val
        return roman_str

    # Total characters allowed on the paper
    paper_capacity = 10000

    # The set of characters Caesar can use in his message
    # (A-Z and space)
    possible_characters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ "
    
    # Initialize minimum cost to a very high value
    min_cost_per_char = float('inf')

    # Iterate through each possible character to find its encryption cost
    # We assume A=1, B=2, ..., Z=26, space=27
    for i in range(len(possible_characters)):
        # The numeric value for the character
        numeric_value = i + 1
        
        # Get the Roman numeral string for this number
        roman_numeral = int_to_roman(numeric_value)
        
        # The cost is the length of the resulting Roman numeral string
        cost = len(roman_numeral)
        
        # Update the minimum cost if the current character is cheaper
        if cost < min_cost_per_char:
            min_cost_per_char = cost

    # The maximum message length is the total capacity divided by the minimum cost
    max_message_length = paper_capacity // min_cost_per_char

    # As requested, print the equation showing how the result is derived.
    # Each number in the final equation is outputted.
    print(f"{paper_capacity} / {min_cost_per_char} = {max_message_length}")

if __name__ == "__main__":
    main()