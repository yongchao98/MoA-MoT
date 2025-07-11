import math

# Define the constants based on the problem description.
total_symbols_in_ching = 9999
num_instructions = 108
num_unique_chars = 108

# Define the size of each part of an instruction in decimal digits (D).
size_of_opcode = 1
size_of_register = 1
size_of_operand = 4

# Step 1: Calculate the size of a single instruction and the total code size.
instruction_size = size_of_opcode + size_of_register + size_of_operand
code_size = num_instructions * instruction_size

# Step 2: Calculate the size of the data portion available for the message.
data_size = total_symbols_in_ching - code_size

# Step 3: Determine the number of decimal digits required to represent one character.
# We need to find the smallest integer 'd' such that 10^d >= num_unique_chars.
# This is equivalent to ceiling(log10(num_unique_chars)).
digits_per_char = math.ceil(math.log10(num_unique_chars))

# Step 4: Calculate K, the highest number of characters that can be decoded.
# This is an integer division since we can't have a fraction of a character.
K = data_size // digits_per_char

# Print the final calculation step-by-step as requested.
print(f"The total space in the Ching is {total_symbols_in_ching} symbols.")
print(f"The program code occupies {num_instructions} instructions * {instruction_size} symbols/instruction = {code_size} symbols.")
print(f"The remaining space for data is {total_symbols_in_ching} - {code_size} = {data_size} symbols.")
print(f"Each of the {num_unique_chars} characters requires {digits_per_char} symbols for representation.")
print("\nThe final equation for the highest number of characters (K) is:")
print(f"K = Data Size / Symbols Per Character")
print(f"K = {data_size} / {digits_per_char} = {K}")

print("<<<3117>>>")