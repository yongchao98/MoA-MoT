def solve_caesar_cipher_length():
    """
    Calculates the length of the longest possible message for Caesar's new encryption.
    """
    print("Thinking Process:")
    print("1. The encryption maps each character (A-Z, space) to a number (1-26, 27).")
    print("2. Each number is converted to a Roman numeral string.")
    print("3. The final encrypted text is the concatenation of these Roman numerals.")
    print("4. The total length of the encrypted text cannot exceed 10000 characters.")
    print("5. To maximize the original message length, we must use characters with the shortest Roman numeral representation (minimum cost).\n")

    def int_to_roman(num):
        """Converts an integer to a Roman numeral string."""
        val_map = [
            (1000, "M"), (900, "CM"), (500, "D"), (400, "CD"),
            (100, "C"), (90, "XC"), (50, "L"), (40, "XL"),
            (10, "X"), (9, "IX"), (5, "V"), (4, "IV"),
            (1, "I")
        ]
        roman_num = []
        for val, sym in val_map:
            while num >= val:
                roman_num.append(sym)
                num -= val
        return "".join(roman_num)

    min_cost = float('inf')
    cheapest_chars = []
    
    # Character set includes A-Z and space
    # A=1, B=2, ..., Z=26, space=27
    num_of_chars = 27 

    print("Calculating the encryption cost (length) for each character:")
    for i in range(1, num_of_chars + 1):
        if i <= 26:
            char = chr(ord('A') + i - 1)
        else:
            char = 'space'
        
        roman_representation = int_to_roman(i)
        cost = len(roman_representation)
        
        # Uncomment the line below to see the cost for each character
        # print(f"  - Character '{char}' ({i}) -> '{roman_representation}' (Cost: {cost})")

        if cost < min_cost:
            min_cost = cost
            cheapest_chars = [char]
        elif cost == min_cost:
            cheapest_chars.append(char)

    paper_limit = 10000
    
    # The longest message is composed of characters with the minimum encryption cost.
    # max_length = total_limit / cost_per_character
    max_message_length = paper_limit // min_cost

    print(f"\nThe characters with the minimum encryption cost of {min_cost} are: {', '.join(cheapest_chars)}")
    print("\nTo write the longest message, Caesar should only use these characters.")
    
    print("\nFinal Calculation:")
    print(f"Paper character limit: {paper_limit}")
    print(f"Minimum characters needed for one letter (cost): {min_cost}")
    print(f"The final equation is: {paper_limit} / {min_cost} = {max_message_length}")
    
    print(f"\nThe length of his longest message is {max_message_length}.")

solve_caesar_cipher_length()
<<<10000>>>