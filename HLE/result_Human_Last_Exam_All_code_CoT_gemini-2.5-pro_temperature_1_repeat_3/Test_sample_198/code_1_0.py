import math

# Step 1: Define the constants from the problem description.
total_symbols_in_ching = 9999
num_instructions = 108
num_unique_characters = 108

# Size of each component of an instruction in decimal digits (D)
opcode_size = 1
register_size = 1
operand_size = 4

# Step 2: Calculate the total size of an instruction and the total size of the code.
digits_per_instruction = opcode_size + register_size + operand_size
total_code_size = num_instructions * digits_per_instruction

# Step 3: Calculate the total size of the data section.
total_data_size = total_symbols_in_ching - total_code_size

# Step 4: Determine the minimum number of digits to represent one of the 108 characters.
# We need to find n such that 10^(n-1) < 108 <= 10^n.
# 10^1 = 10 (< 108)
# 10^2 = 100 (< 108)
# 10^3 = 1000 (>= 108)
# So, we need 3 digits per character.
digits_per_character = math.ceil(math.log10(num_unique_characters))

# Step 5: Calculate the highest number of characters (K) that can be decoded.
# This is the total data size divided by the number of digits needed for each character.
# We use integer division because we can only decode a whole number of characters.
K = total_data_size // digits_per_character

# Final output: Print the equation with all the numbers used.
print(f"The highest number of characters K can be calculated as:")
print(f"K = (Total Symbols - Code Size) / Digits per Character")
print(f"K = ({total_symbols_in_ching} - {num_instructions} * {digits_per_instruction}) // {digits_per_character}")
print(f"K = ({total_symbols_in_ching} - {total_code_size}) // {digits_per_character}")
print(f"K = {total_data_size} // {digits_per_character}")
print(f"K = {K}")