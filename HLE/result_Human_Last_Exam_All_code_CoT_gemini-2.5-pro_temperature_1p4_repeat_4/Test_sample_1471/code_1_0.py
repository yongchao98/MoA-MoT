def roman_to_int(s):
    """Converts a Roman numeral string to an integer."""
    roman_map = {'I': 1, 'V': 5, 'X': 10, 'L': 50, 'C': 100, 'D': 500, 'M': 1000}
    result = 0
    prev_value = 0
    # Iterate through the string from right to left for easier conversion logic
    for char in reversed(s):
        current_value = roman_map[char]
        # If the current value is less than the previous one, it's a subtractive case (e.g., IV, IX)
        if current_value < prev_value:
            result -= current_value
        else:
            result += current_value
        prev_value = current_value
    return result

# The Roman numerals from the original equation
roman_numeral_1 = 'VIII'
roman_numeral_2 = 'CCC'
roman_numeral_3 = 'CCCI'

# Convert each Roman numeral to its modern integer equivalent
num1 = roman_to_int(roman_numeral_1)
num2 = roman_to_int(roman_numeral_2)
num3 = roman_to_int(roman_numeral_3)

# Perform the calculation
solution = num1 / num2 - num3

# Print the equation with modern numbers and the final result
print(f"The equation with modern numbers is:")
print(f"{num1} / {num2} - {num3} = {solution}")