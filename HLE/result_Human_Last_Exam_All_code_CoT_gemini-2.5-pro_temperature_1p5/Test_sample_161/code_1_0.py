def solve_caesar_cipher():
    """
    Calculates the length of the longest possible message Caesar can write.
    """
    # The total number of characters the paper can hold for the encrypted message.
    paper_limit = 10000

    def int_to_roman(num):
        """Converts an integer to its Roman numeral representation."""
        val = [
            1000, 900, 500, 400,
            100, 90, 50, 40,
            10, 9, 5, 4,
            1
        ]
        syb = [
            "M", "CM", "D", "CD",
            "C", "XC", "L", "XL",
            "X", "IX", "V", "IV",
            "I"
        ]
        roman_num = ''
        i = 0
        while num > 0:
            for _ in range(num // val[i]):
                roman_num += syb[i]
                num -= val[i]
            i += 1
        return roman_num

    # The problem involves letters A-Z and space. We map them to numbers 1-27.
    # A=1, B=2, ..., Z=26, space=27.
    character_numbers = range(1, 28)

    # Calculate the encryption cost (length of Roman numeral) for each possible character.
    min_cost = float('inf')
    for num in character_numbers:
        roman_representation = int_to_roman(num)
        cost = len(roman_representation)
        if cost < min_cost:
            min_cost = cost

    # The maximum message length is the paper limit divided by the minimum cost.
    # We use integer division as we can only have whole characters.
    max_length = paper_limit // min_cost

    # Print the explanation and the final equation.
    print("To maximize the message length, we must use characters with the shortest Roman numeral representation.")
    print(f"The minimum cost to encrypt a single character is {min_cost}.")
    print("The total character limit on the paper for the encrypted message is 10000.")
    print("\nThe final equation to calculate the maximum message length is:")
    print(f"{paper_limit} / {min_cost} = {max_length}")
    print(f"\nTherefore, the longest possible message is {max_length} characters.")


solve_caesar_cipher()
<<<10000>>>