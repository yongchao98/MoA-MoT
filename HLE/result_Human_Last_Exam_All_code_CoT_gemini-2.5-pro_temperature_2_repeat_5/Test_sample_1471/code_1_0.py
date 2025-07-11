def roman_to_int(s):
    """Converts a Roman numeral string to an integer."""
    roman_map = {'I': 1, 'V': 5, 'X': 10, 'L': 50, 'C': 100, 'D': 500, 'M': 1000}
    
    # Handle subtractive notation by replacing with additive equivalents
    s = s.replace("IV", "IIII").replace("IX", "VIIII")
    s = s.replace("XL", "XXXX").replace("XC", "LXXXX")
    s = s.replace("CD", "CCCC").replace("CM", "DCCCC")
    
    num = 0
    for char in s:
        num += roman_map[char]
    return num

# The Roman numerals from the equation
roman1 = "VIII"
roman2 = "CCC"
roman3 = "CCLIII"

# Convert them to integers
num1 = roman_to_int(roman1)
num2 = roman_to_int(roman2)
num3 = roman_to_int(roman3)

# Perform the calculation
result = (num1 / num2) - num3

# Print the equation in modern numbers and the result
print(f"The equation {roman1} / {roman2} - {roman3} translates to:")
print(f"{num1} / {num2} - {num3} = {result}")
