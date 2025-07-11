import math

# Step 1: Define constants and calculate code size.
num_instructions = 108
instruction_size_d = 1 + 1 + 4  # opcode (1) + register (1) + operand (4)
code_size_d = num_instructions * instruction_size_d

# Step 2: Determine the boundaries of the available data space.
# Code is loaded at address 0.
data_start_address = code_size_d
# The last usable memory address is 9997, as 9998 and 9999 are for I/O.
data_end_address = 9997

# Calculate the total size of the memory space available for the message data.
data_space_d = data_end_address - data_start_address + 1

# Step 3: Define the size of the data block for a single character.
data_per_char_d = 4

# Step 4: Calculate the highest number of characters, K.
# We use integer division (or floor) since we can't have a fraction of a character.
K = data_space_d // data_per_char_d

# Output the final calculation, showing all the numbers involved in the equation.
print(f"To find K, we calculate the available data space and divide it by the data size per character.")
print(f"Code size = {num_instructions} instructions * {instruction_size_d} digits/instruction = {code_size_d} digits.")
print(f"Available data space = (Last usable address - Start data address + 1) = ({data_end_address} - {data_start_address} + 1) = {data_space_d} digits.")
print(f"Data per character = {data_per_char_d} digits.")
print("\nFinal Equation:")
print(f"K = floor(Available data space / Data per character)")
print(f"K = floor({data_space_d} / {data_per_char_d}) = {K}")
<<<2337>>>