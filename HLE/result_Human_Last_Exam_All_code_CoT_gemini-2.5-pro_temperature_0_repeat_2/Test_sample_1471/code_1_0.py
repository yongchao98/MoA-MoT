def roman_to_int(s):
    """Converts a Roman numeral string to an integer."""
    rom_val = {'I': 1, 'V': 5, 'X': 10, 'L': 50, 'C': 100, 'D': 500, 'M': 1000}
    int_val = 0
    i = 0
    while i < len(s):
        # Check for subtractive notation (e.g., IV, IX, XL, XC, CD, CM)
        if i + 1 < len(s) and rom_val[s[i]] < rom_val[s[i+1]]:
            int_val += rom_val[s[i+1]] - rom_val[s[i]]
            i += 2
        else:
            int_val += rom_val[s[i]]
            i += 1
    return int_val

# Roman numerals from the equation
r1 = "XL"
r2 = "CCC"
r3 = "CCMCMXCIX"

# Convert Roman numerals to integers
n1 = roman_to_int(r1)
n2 = roman_to_int(r2)
n3 = roman_to_int(r3)

# Perform the calculation
result = n1 / n2 - n3

# Print the equation with modern numbers
print(f"The equation is: {n1} / {n2} - {n3} = ?")

# Print the final answer
print(f"The solution is: {result}")

# The final answer in the required format
print(f"<<<{result}>>>")