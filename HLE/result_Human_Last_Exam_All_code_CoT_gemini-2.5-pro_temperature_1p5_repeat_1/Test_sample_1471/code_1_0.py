def roman_to_int(s):
    """
    Converts a Roman numeral string to an integer.
    """
    roman_map = {'I': 1, 'V': 5, 'X': 10, 'L': 50, 'C': 100, 'D': 500, 'M': 1000}
    result = 0
    for i in range(len(s)):
        # If the next numeral is larger, this is a subtractive case (e.g., IV, IX)
        if i + 1 < len(s) and roman_map[s[i]] < roman_map[s[i+1]]:
            result -= roman_map[s[i]]
        else:
            result += roman_map[s[i]]
    return result

# The Roman numerals from the equation
roman_str1 = "VIII"
roman_str2 = "CCC"
roman_str3 = "CCCXI"

# Convert Roman numerals to modern numbers
num1 = roman_to_int(roman_str1)
num2 = roman_to_int(roman_str2)
num3 = roman_to_int(roman_str3)

# Perform the calculation
# The equation is: (VIII / CCC) - CCCXI
solution = (num1 / num2) - num3

# Print the equation with modern numbers and the final answer
print(f"{num1} / {num2} - {num3} = {solution}")
<<<_solution_>>>