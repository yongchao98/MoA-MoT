import math

def solve_caesar_riddle():
    """
    This function calculates the length of the longest possible message Caesar can write.
    """
    
    # The paper's total character capacity for the encrypted message.
    paper_limit = 10000

    # This function converts an integer to its Roman numeral representation.
    def int_to_roman(num):
        if not isinstance(num, int) or not 0 < num < 4000:
            return "" # Roman numerals are typically for positive integers.
        
        val_map = [
            (1000, "M"), (900, "CM"), (500, "D"), (400, "CD"),
            (100, "C"), (90, "XC"), (50, "L"), (40, "XL"),
            (10, "X"), (9, "IX"), (5, "V"), (4, "IV"),
            (1, "I")
        ]
        roman_num = ""
        for val, symbol in val_map:
            # Add the Roman numeral symbol for the largest value less than or equal to the number.
            while num >= val:
                roman_num += symbol
                num -= val
        return roman_num

    # To find the longest message, we must use characters that have the shortest
    # Roman numeral representation (i.e., the lowest "cost").
    # Let's find this minimum cost by checking the letters A-Z (mapped to numbers 1-26).
    
    min_cost_per_char = float('inf')
    
    # Iterate through the numerical values for A (1) to Z (26).
    for i in range(1, 27):
        roman_numeral = int_to_roman(i)
        cost = len(roman_numeral)
        if cost < min_cost_per_char:
            min_cost_per_char = cost
    
    # Note: The 'space' character is not defined, but its encrypted length cannot be less than the minimum we find.
    # The minimum cost will be 1, for characters like 'A' (I), 'E' (V), or 'J' (X).

    # The maximum length of the plaintext is the total paper limit divided by the minimum cost per character.
    max_message_length = math.floor(paper_limit / min_cost_per_char)

    print(f"The total character limit for the encrypted message is {paper_limit}.")
    print(f"To create the longest message, we must use characters that encrypt to the shortest Roman numerals.")
    print(f"The minimum number of Roman numerals required to represent one character (e.g., 'A' -> 'I') is {min_cost_per_char}.")
    print("\nThe equation to find the maximum message length is:")
    print(f"{paper_limit} (total limit) / {min_cost_per_char} (cost per character) = {max_message_length}")
    print(f"\nThus, the length of his longest message is {max_message_length} characters.")

solve_caesar_riddle()
<<<10000>>>