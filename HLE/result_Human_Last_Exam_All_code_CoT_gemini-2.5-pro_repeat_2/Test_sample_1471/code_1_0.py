def roman_to_int(s):
    """Converts a Roman numeral string to an integer."""
    roman_map = {'I': 1, 'V': 5, 'X': 10, 'L': 50, 'C': 100, 'D': 500, 'M': 1000}
    result = 0
    i = 0
    while i < len(s):
        # Check for subtractive notation (e.g., IV, IX)
        if i + 1 < len(s) and roman_map[s[i]] < roman_map[s[i+1]]:
            result += roman_map[s[i+1]] - roman_map[s[i]]
            i += 2
        else:
            result += roman_map[s[i]]
            i += 1
    return result

# Roman numerals from the equation
roman1 = 'IV'
roman2 = 'XXX'
roman3 = 'XXXI'

# Convert Roman numerals to modern numbers
num1 = roman_to_int(roman1)
num2 = roman_to_int(roman2)
num3 = roman_to_int(roman3)

# Perform the calculation
solution = num1 / num2 - num3

# Print the final equation with the numbers and the result
print(f"{num1} / {num2} - {num3} = {solution}")

# Also print the final answer in the requested format
# print(f"<<<{solution}>>>")