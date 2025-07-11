import math

# Step 1: Define the constants from the problem description.
ching_total_digits = 9999
num_instructions = 108
opcode_digits = 1
register_digits = 1
operand_digits = 4

# The fundamental unit of data processing is determined by the register size.
data_chunk_size = 4 # Registers are 4D, so data is processed in 4-digit chunks.

# Step 2: Calculate the total size of the program code in digits.
# Each instruction has a fixed size.
instruction_size = opcode_digits + register_digits + operand_digits
# The total program size is the number of instructions multiplied by the size of each one.
program_size = num_instructions * instruction_size

# Step 3: Calculate the total size of the encrypted data.
# The data is the part of the Ching that is not the program.
encrypted_data_size = ching_total_digits - program_size

# Step 4: Calculate K, the highest number of characters that can be decoded.
# This is the total data size divided by the amount of data needed per character.
# We use integer division (floor) since we can't decode a fraction of a character.
K = encrypted_data_size // data_chunk_size

# Final Output: Print the reasoning and the final equation with all values.
print("This script calculates K, the highest number of characters that can be decoded.")
print("-" * 50)
print(f"1. Calculate total program size:")
print(f"   ({opcode_digits} + {register_digits} + {operand_digits}) digits/instruction * {num_instructions} instructions = {program_size} digits")
print("-" * 50)
print(f"2. Calculate total encrypted data size:")
print(f"   {ching_total_digits} total digits - {program_size} program digits = {encrypted_data_size} data digits")
print("-" * 50)
print(f"3. Calculate K based on processing data in {data_chunk_size}-digit chunks:")
print(f"   K = floor({encrypted_data_size} data digits / {data_chunk_size} digits/character)")
print(f"   K = {K}")
print("-" * 50)
print("The final equation for K is:")
print(f"K = ({ching_total_digits} - ({num_instructions} * ({opcode_digits} + {register_digits} + {operand_digits}))) // {data_chunk_size}")
print(f"K = {K}")
