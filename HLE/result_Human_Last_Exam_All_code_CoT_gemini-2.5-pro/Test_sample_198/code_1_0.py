import math

# Step 1: Define constants and calculate the program code size.
# The problem describes the memory, instructions, and characters.
total_memory_D = 10000  # Total memory in Decimal Digits (D).
num_instructions = 108  # The program has 108 instructions.
# Each instruction is [opcode][register][operand] = 1D + 1D + 4D.
instruction_size_D = 1 + 1 + 4
code_size_D = num_instructions * instruction_size_D

# Step 2: Calculate the available data space.
# Two memory locations (9999 and 9998) are for I/O.
io_buffer_size_D = 2
data_space_D = total_memory_D - code_size_D - io_buffer_size_D

# Step 3: Determine the number of digits needed per character.
# There are 108 unique characters from the Shu.
num_unique_characters = 108
# To represent 108 unique items in a decimal system (base 10), we need
# to find the smallest integer 'd' where 10^d >= 108.
# This is equivalent to ceil(log10(108)).
digits_per_char = math.ceil(math.log10(num_unique_characters))

# Step 4: Calculate the highest number of characters, K.
# This is the total data space divided by the space needed for one character.
# We use integer division because we can only decode whole characters.
K = data_space_D // digits_per_char

# Print the step-by-step calculation.
print(f"1. Calculate Code Size:")
print(f"   - Number of instructions: {num_instructions}")
print(f"   - Size of one instruction: {instruction_size_D}D")
print(f"   - Total code size = {num_instructions} * {instruction_size_D} = {code_size_D}D\n")

print(f"2. Calculate Available Data Space:")
print(f"   - Total memory: {total_memory_D}D")
print(f"   - Reserved for I/O: {io_buffer_size_D}D")
print(f"   - Data space = {total_memory_D} - {code_size_D} - {io_buffer_size_D} = {data_space_D}D\n")

print(f"3. Determine Digits per Character:")
print(f"   - Number of unique characters: {num_unique_characters}")
print(f"   - To represent {num_unique_characters} characters, we need ceil(log10({num_unique_characters})) = {digits_per_char} decimal digits per character.\n")

print(f"4. Calculate K (Highest Number of Characters):")
print(f"   - K = floor(Available Data Space / Digits per Character)")
print(f"   - K = floor({data_space_D} / {digits_per_char})")
print(f"   - K = {K}\n")

print(f"Final Equation:")
print(f"K = ({total_memory_D} - ({num_instructions} * {instruction_size_D}) - {io_buffer_size_D}) // {digits_per_char} = {K}")