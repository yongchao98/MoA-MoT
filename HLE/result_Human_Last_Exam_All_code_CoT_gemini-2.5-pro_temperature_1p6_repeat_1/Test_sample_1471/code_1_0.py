def roman_to_int(s):
    """Converts a Roman numeral string to an integer."""
    roman_map = {'I': 1, 'V': 5, 'X': 10, 'L': 50, 'C': 100, 'D': 500, 'M': 1000}
    result = 0
    for i in range(len(s)):
        # Check for subtractive notation (e.g., IV, IX, XL)
        if i > 0 and roman_map[s[i]] > roman_map[s[i-1]]:
            result += roman_map[s[i]] - 2 * roman_map[s[i-1]]
        else:
            result += roman_map[s[i]]
    return result

# The Roman numerals from the equation
roman_1 = "VIII"
roman_2 = "CCC"
roman_3 = "CCXII"

# Convert Roman numerals to modern numbers
num_1 = roman_to_int(roman_1)
num_2 = roman_to_int(roman_2)
num_3 = roman_to_int(roman_3)

# Perform the calculation
solution = num_1 / num_2 - num_3

# Print the equation with modern numbers and the solution
print(f"{num_1} / {num_2} - {num_3} = {solution}")