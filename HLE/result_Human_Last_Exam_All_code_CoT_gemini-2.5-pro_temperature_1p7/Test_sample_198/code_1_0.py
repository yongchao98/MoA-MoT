import math

# Define the constants from the problem description
total_symbols = 9999
num_instructions = 108
num_unique_characters = 108

# Calculate the size of a single instruction in decimal digits (D)
# Format: [opcode (1D)] + [register (1D)] + [operand (4D)]
opcode_size = 1
register_size = 1
operand_size = 4
instruction_size_d = opcode_size + register_size + operand_size

# Calculate the total size of the program code in digits
code_size_d = num_instructions * instruction_size_d

# Calculate the total size of the data section in digits
data_size_d = total_symbols - code_size_d

# Calculate the number of digits required to represent a single Chinese character.
# We need 'n' digits such that 10^n >= 108.
# n = ceil(log10(108))
digits_per_character = math.ceil(math.log10(num_unique_characters))

# Calculate the highest number of characters (K) that can be decoded
# This is the integer division of the total data size by the size per character.
K = data_size_d // digits_per_character

# Print the step-by-step calculation
print("Step 1: Calculate the total size of the program code.")
print(f"Size of one instruction = {opcode_size} (opcode) + {register_size} (register) + {operand_size} (operand) = {instruction_size_d} digits.")
print(f"Total code size = {num_instructions} instructions * {instruction_size_d} digits/instruction = {code_size_d} digits.")
print("-" * 30)

print("Step 2: Calculate the available space for data.")
print(f"Total data space = {total_symbols} (total digits) - {code_size_d} (code digits) = {data_size_d} digits.")
print("-" * 30)

print("Step 3: Calculate the digits needed per character.")
print(f"To represent {num_unique_characters} unique items, we need ceil(log10({num_unique_characters})) = {digits_per_character} decimal digits per item.")
print("-" * 30)

print("Step 4: Calculate the highest number of characters (K).")
print(f"K = Total Data Space / Digits per Character")
print(f"K = {data_size_d} / {digits_per_character} = {K}")
print(f"\nThe highest number of characters that can be decoded is {K}.")

<<<3117>>>