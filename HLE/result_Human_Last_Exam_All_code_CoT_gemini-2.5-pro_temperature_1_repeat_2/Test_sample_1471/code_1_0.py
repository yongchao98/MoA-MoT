def roman_to_int(s):
    """
    Converts a Roman numeral string to an integer.
    """
    roman_map = {'I': 1, 'V': 5, 'X': 10, 'L': 50, 'C': 100, 'D': 500, 'M': 1000}
    result = 0
    i = 0
    while i < len(s):
        # Check for subtractive notation (e.g., IV, IX, XL)
        if i + 1 < len(s) and roman_map[s[i]] < roman_map[s[i+1]]:
            result += roman_map[s[i+1]] - roman_map[s[i]]
            i += 2
        else:
            result += roman_map[s[i]]
            i += 1
    return result

# The Roman numerals from the equation
roman_num1 = "VIII"
roman_num2 = "CCC"
roman_num3 = "CCXCVIII"

# Convert Roman numerals to modern integers
num1 = roman_to_int(roman_num1)
num2 = roman_to_int(roman_num2)
num3 = roman_to_int(roman_num3)

# Perform the calculation
# Equation: num1 / num2 - num3
solution = (num1 / num2) - num3

# Print the equation in modern numbers with the solution
print(f"{num1} / {num2} - {num3} = {solution}")
