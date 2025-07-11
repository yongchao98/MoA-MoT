def int_to_roman(num):
    """Converts an integer to a Roman numeral string."""
    # Standard values and symbols for Roman numerals
    val = [1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1]
    syb = ["M", "CM", "D", "CD", "C", "XC", "L", "XL", "X", "IX", "V", "IV", "I"]
    roman_num = ""
    i = 0
    while num > 0:
        for _ in range(num // val[i]):
            roman_num += syb[i]
            num -= val[i]
        i += 1
    return roman_num

def solve_longest_message():
    """
    Calculates the length of the longest possible message based on the described
    Roman numeral encryption.
    """
    # Step 1 & 2: Define character values and find the minimum Roman numeral length.
    # We consider all 27 possible characters (A-Z and space) mapping to values 1-27.
    char_values = range(1, 28)

    min_roman_len = float('inf')
    cheapest_char_value = 0

    for value in char_values:
        roman_len = len(int_to_roman(value))
        if roman_len < min_roman_len:
            min_roman_len = roman_len
            cheapest_char_value = value

    # Step 3: Calculate the minimum cost per character.
    # To avoid ambiguity, each encoded character requires a separator of length 1.
    separator_len = 1
    min_cost_per_char = min_roman_len + separator_len
    
    # Step 4 & 5: Calculate the maximum message length.
    paper_capacity = 10000
    max_length = paper_capacity // min_cost_per_char

    print("--- Analysis ---")
    print(f"1. The shortest possible Roman numeral for a character (A-Z, space) corresponds to the number {cheapest_char_value}, which is '{int_to_roman(cheapest_char_value)}'.")
    print(f"2. The length of this shortest Roman numeral is {min_roman_len}.")
    print(f"3. To avoid ambiguity, a separator of length {separator_len} is added after each character's code.")
    print(f"4. The minimum cost to write one character is the sum of these lengths.")
    print(f"   Minimum Cost = {min_roman_len} (Roman Numeral) + {separator_len} (Separator) = {min_cost_per_char}")
    print("\n--- Final Calculation ---")
    print("The total capacity of the paper is 10000 characters.")
    print("The final equation for the longest message length is: Total Capacity / Minimum Cost")
    print(f"Longest Message Length = {paper_capacity} / {min_cost_per_char}")
    print(f"Result: {max_length}")

if __name__ == '__main__':
    solve_longest_message()
