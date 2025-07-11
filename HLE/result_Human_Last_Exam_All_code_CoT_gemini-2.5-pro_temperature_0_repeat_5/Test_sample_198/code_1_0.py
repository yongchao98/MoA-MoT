import math

# Step 1: Define the constants from the problem description.
# Total symbols (decimal digits) in the Ching.
total_symbols = 9999
# Number of instructions in the program.
num_instructions = 108
# Size of a single instruction: 1D for opcode, 1D for register, 4D for operand.
instruction_size = 1 + 1 + 4
# Number of unique Chinese characters the program can output.
num_unique_chars = 108

# Step 2: Calculate the size of the program code and the remaining data.
# The code section is the number of instructions multiplied by the size of each.
code_size = num_instructions * instruction_size
# The data section is the total symbols minus the code size.
data_size = total_symbols - code_size

# Step 3: Determine the minimum number of digits to represent one character.
# We need to find the smallest integer S such that 10^S >= num_unique_chars.
# This is equivalent to ceil(log10(num_unique_chars)).
digits_per_char = math.ceil(math.log10(num_unique_chars))

# Step 4: Calculate K, the maximum number of characters.
# This is the integer division of the total data size by the digits needed per character.
K = data_size // digits_per_char

# Final Step: Print the equation showing how K is calculated, including all the numbers.
print(f"K = ({total_symbols} - ({num_instructions} * {instruction_size})) / {digits_per_char} = {K}")

# The final answer for K
# print(f"<<<{K}>>>")