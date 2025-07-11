def roman_to_int(s: str) -> int:
    """Converts a Roman numeral string to an integer."""
    roman_map = {'I': 1, 'V': 5, 'X': 10, 'L': 50, 'C': 100, 'D': 500, 'M': 1000}
    result = 0
    i = 0
    while i < len(s):
        # Check for the subtractive case (e.g., IX for 9 or IV for 4)
        if i + 1 < len(s) and roman_map[s[i]] < roman_map[s[i+1]]:
            result += roman_map[s[i+1]] - roman_map[s[i]]
            i += 2
        # Otherwise, handle the additive case (e.g., VI for 6 or XI for 11)
        else:
            result += roman_map[s[i]]
            i += 1
    return result

# The Roman numerals from the equation "VIII/LXXX - XXI/IX = ?"
roman_num_1 = "VIII"
roman_num_2 = "LXXX"
roman_num_3 = "XXI"
roman_num_4 = "IX"

# Convert the Roman numerals to integers
num_1 = roman_to_int(roman_num_1)
num_2 = roman_to_int(roman_num_2)
num_3 = roman_to_int(roman_num_3)
num_4 = roman_to_int(roman_num_4)

# Calculate the result of the equation
solution = (num_1 / num_2) - (num_3 / num_4)

# Print the full equation with modern numbers and the final answer
print(f"The equation is: {num_1} / {num_2} - {num_3} / {num_4} = {solution}")