import sys

def int_to_roman(num):
    """Converts an integer to a Roman numeral string."""
    val_map = [
        (1000, "M"), (900, "CM"), (500, "D"), (400, "CD"),
        (100, "C"), (90, "XC"), (50, "L"), (40, "XL"),
        (10, "X"), (9, "IX"), (5, "V"), (4, "IV"),
        (1, "I")
    ]
    roman_num = []
    for val, numeral in val_map:
        count = num // val
        roman_num.append(numeral * count)
        num %= val
    return "".join(roman_num)

def solve_longest_message():
    """
    Calculates the length of the longest possible message given the encryption rules.
    """
    # The paper's total character capacity.
    total_capacity = 10000

    # We need to find the minimum cost to write a single character. A character's
    # cost is the length of its Roman numeral representation. We'll map A=1, B=2, ..., Z=26.
    min_cost = sys.maxsize
    
    # Calculate the cost for each letter from A (mapped to 1) to Z (mapped to 26).
    for i in range(1, 27):
        roman_representation = int_to_roman(i)
        cost = len(roman_representation)
        if cost < min_cost:
            min_cost = cost
            
    # The cost of the 'space' character is not specified. However, any Roman numeral
    # representation must have a length of at least 1. The letters 'A' (1 -> "I"), 
    # 'E' (5 -> "V"), and 'J' (10 -> "X") all have a cost of 1.
    # Therefore, the minimum possible cost for any character is 1.

    # To write the longest message, Caesar should use only characters with this minimum cost.
    # The maximum length of the message is the total capacity divided by the minimum cost.
    max_length = total_capacity // min_cost

    print("To find the longest message, we must use the characters with the shortest encrypted length (cost).")
    print("The cost of a character is the length of its Roman numeral representation.")
    print(f"\nTotal capacity of the paper: {total_capacity} characters.")
    print(f"Minimum cost per character: {min_cost} (e.g., for 'A' -> 'I')")
    
    print("\nThe final equation to find the maximum message length is:")
    # Outputting each number in the final equation as requested.
    print(f"{total_capacity} / {min_cost} = {max_length}")
    
    print(f"\nTherefore, the length of the longest possible message is {max_length}.")

solve_longest_message()
<<<10000>>>