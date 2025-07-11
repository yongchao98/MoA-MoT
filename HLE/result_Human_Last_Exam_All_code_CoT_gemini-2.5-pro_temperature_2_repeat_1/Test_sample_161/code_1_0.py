import math

def integer_to_roman(num):
    """Converts a positive integer to its Roman numeral representation."""
    val_map = [
        (1000, "M"), (900, "CM"), (500, "D"), (400, "CD"),
        (100, "C"), (90, "XC"), (50, "L"), (40, "XL"),
        (10, "X"), (9, "IX"), (5, "V"), (4, "IV"),
        (1, "I")
    ]
    roman_num = []
    for val, numeral in val_map:
        while num >= val:
            roman_num.append(numeral)
            num -= val
    return "".join(roman_num)

# There are 27 characters in the alphabet (A-Z and space).
# We find the shortest Roman numeral representations by checking numbers 1 to 27.
num_chars_in_alphabet = 27
roman_lengths = []
for i in range(1, num_chars_in_alphabet + 1):
    roman_lengths.append(len(integer_to_roman(i)))

# The minimum length of a Roman numeral for numbers 1-27
min_roman_length = min(roman_lengths)

# The total capacity of the paper
paper_capacity = 10000

# The length of the delimiter is 1 character.
delimiter_length = 1

# Let n be the length of the longest possible message.
# The total encrypted length is: n * min_roman_length + (n - 1) * delimiter_length
# This must be <= paper_capacity
# n * min_roman_length + n - 1 <= 10000
# n * (min_roman_length + 1) - 1 <= 10000
# n * (min_roman_length + 1) <= 10001
# n <= 10001 / (min_roman_length + 1)
# We need the floor of this value for the maximum integer length.

# We print the steps of the calculation as requested.
print(f"To create a decodable message, each Roman numeral must be separated by a delimiter (e.g., a comma or space).")
print(f"Let the length of the message be n.")
print(f"The number of characters to encode is: {num_chars_in_alphabet}")
print(f"The minimum possible length of a Roman numeral for numbers 1-{num_chars_in_alphabet} is: {min_roman_length}")
print(f"The length of the delimiter is: {delimiter_length}")
print(f"The paper has a capacity of: {paper_capacity} characters.")
print(f"The equation for the total encrypted length is: n * {min_roman_length} + (n - 1) * {delimiter_length} <= {paper_capacity}")
print(f"This simplifies to: {min_roman_length + delimiter_length}*n - {delimiter_length} <= {paper_capacity}")
numerator = paper_capacity + delimiter_length
denominator = min_roman_length + delimiter_length
print(f"Solving for n: n <= ({paper_capacity} + {delimiter_length}) / ({min_roman_length} + {delimiter_length})")
print(f"So, n <= {numerator} / {denominator}")

max_length = math.floor(numerator / denominator)

print(f"\nThe maximum integer length of the message is the floor of the result.")
print(f"Maximum message length = {max_length}")
