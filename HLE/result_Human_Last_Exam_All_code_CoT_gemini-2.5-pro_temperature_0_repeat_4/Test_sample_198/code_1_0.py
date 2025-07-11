import math

# This script calculates K, the highest number of characters that can be decoded.

# Define the constants based on the information provided in the problem.
total_symbols_in_ching = 9999
num_instructions = 108
num_unique_characters = 108

# Each instruction consists of [opcode][register][operand], which is 1 + 1 + 4 = 6 decimal digits.
digits_per_instruction = 6

# 1. Calculate the total size of the program code in decimal digits.
code_size_in_digits = num_instructions * digits_per_instruction

# 2. Calculate the total size of the data section.
data_size_in_digits = total_symbols_in_ching - code_size_in_digits

# 3. Determine the number of decimal digits needed to represent one of the 108 unique characters.
# We need to find the smallest integer 'd' such that 10^d >= 108.
# This can be calculated using the ceiling of the base-10 logarithm.
digits_per_character = math.ceil(math.log10(num_unique_characters))

# 4. Calculate the maximum number of characters (K) that can be encoded in the data section.
# We use integer division because we can only decode whole characters.
K = data_size_in_digits // digits_per_character

# Print the final equation showing how K is derived from the calculated values.
print("The final equation to find K, the highest number of characters, is:")
print(f"K = (Total Symbols - Code Size) / Digits per Character")
print(f"K = ({total_symbols_in_ching} - {num_instructions} * {digits_per_instruction}) / {digits_per_character}")
print(f"K = {data_size_in_digits} / {digits_per_character} = {K}")