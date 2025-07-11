import math

# Define the constants based on the problem description
total_symbols = 9999
num_instructions = 108
num_unique_characters = 108

# Define the size of each part of an instruction in decimal digits (D)
opcode_size = 1
register_size = 1
operand_size = 4

# Calculate the total size of a single instruction
digits_per_instruction = opcode_size + register_size + operand_size

# Calculate the total size of the program's code segment
code_size = num_instructions * digits_per_instruction

# Calculate the total size of the data segment available for the secret message
data_size = total_symbols - code_size

# Calculate the number of digits required to represent a single character.
# This is the smallest integer 'n' where 10^n >= 108.
# We can calculate this using the ceiling of the base-10 logarithm.
digits_per_character = math.ceil(math.log10(num_unique_characters))

# Calculate K, the highest number of characters that can be decoded.
# This is the integer division of the data size by the digits needed per character.
K = data_size // digits_per_character

# Print the final equation showing all the numbers involved in the calculation.
# This equation represents: K = (Total Data) / (Encoding size per character)
print(f"K = ({total_symbols} - {num_instructions} * {digits_per_instruction}) / {digits_per_character} = {K}")