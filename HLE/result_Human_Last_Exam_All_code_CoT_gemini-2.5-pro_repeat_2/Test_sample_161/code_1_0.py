def solve_caesar_cipher():
    """
    Calculates the length of the longest possible message Caesar can write
    under the given constraints.
    """

    def int_to_roman(num):
        """Converts an integer to its Roman numeral representation."""
        # A standard mapping of values to Roman numerals
        val_map = [
            (1000, "M"), (900, "CM"), (500, "D"), (400, "CD"),
            (100, "C"), (90, "XC"), (50, "L"), (40, "XL"),
            (10, "X"), (9, "IX"), (5, "V"), (4, "IV"), (1, "I")
        ]
        roman_numeral = []
        for val, numeral in val_map:
            # Greedily append the largest possible numeral
            count = num // val
            roman_numeral.append(numeral * count)
            num -= val * count
        return "".join(roman_numeral)

    paper_capacity = 10000
    min_cost = float('inf')

    # To maximize the message length, we must find the character
    # with the shortest Roman numeral representation (the minimum cost).
    # We check the costs for letters A (1) through Z (26).
    for i in range(1, 27):
        roman_numeral = int_to_roman(i)
        cost = len(roman_numeral)
        if cost < min_cost:
            min_cost = cost

    # The minimum cost is 1, for characters like 'A' (I), 'E' (V), or 'J' (X).
    # The longest message is the total capacity divided by this minimum cost.
    max_length = paper_capacity // min_cost

    print("To write the longest message, Caesar must use characters that encrypt to the shortest Roman numerals.")
    print(f"The minimum length (cost) of an encrypted character is {min_cost}.")
    print("This is achieved with letters like 'A' (I), 'E' (V), or 'J' (X).")
    print("\nThe maximum message length is calculated by dividing the paper capacity by this minimum cost.")
    print("Final calculation:")
    print(f"{paper_capacity} / {min_cost} = {max_length}")

solve_caesar_cipher()
<<<10000>>>