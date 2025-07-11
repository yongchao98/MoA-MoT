def roman_to_int(s):
    """
    Converts a Roman numeral string to an integer.
    """
    roman_map = {'I': 1, 'V': 5, 'X': 10, 'L': 50, 'C': 100, 'D': 500, 'M': 1000}
    result = 0
    i = 0
    while i < len(s):
        # Get value of the current symbol
        s1 = roman_map[s[i]]

        if (i + 1) < len(s):
            # Get value of the next symbol
            s2 = roman_map[s[i + 1]]

            # Compare symbols
            if s1 >= s2:
                # Value of current symbol is greater or equal to the next symbol
                result = result + s1
                i = i + 1
            else:
                # Value of current symbol is less than the next symbol
                result = result + s2 - s1
                i = i + 2
        else:
            result = result + s1
            i = i + 1
    return result

# The Roman numerals from the equation
roman_a = "VIII"
roman_b = "CCC"
roman_c = "CCXIII"

# Convert Roman numerals to integers
num_a = roman_to_int(roman_a)
num_b = roman_to_int(roman_b)
num_c = roman_to_int(roman_c)

# Calculate the result
solution = (num_a / num_b) - num_c

# Print the equation with modern numbers and the solution
print(f"The equation is: {num_a} / {num_b} - {num_c} = {solution}")
print(f"The original Roman numeral equation was: {roman_a} / {roman_b} - {roman_c}")
print(f"The final answer is: {solution}")
