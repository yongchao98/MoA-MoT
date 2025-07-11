import math

# This script calculates K, the highest number of characters that can be decoded from the Ching.

# --- Step 1: Define constants from the problem description ---
total_symbols = 9999
num_instructions = 108
num_unique_chars = 108

# The size of a single instruction in decimal digits (D)
# Format: [opcode][register][operand] -> 1D + 1D + 4D
instruction_size_d = 1 + 1 + 4

# --- Step 2: Calculate the size of the program code and the remaining data ---
# The program code is loaded first in memory.
total_code_size_d = num_instructions * instruction_size_d

# The rest of the memory is data.
total_data_size_d = total_symbols - total_code_size_d

# --- Step 3: Determine the number of digits needed to represent one character ---
# To represent N unique items using a base-B system, you need ceil(log_B(N)) digits.
# Here, we have 108 characters (N=108) and use the decimal system (B=10).
digits_per_char = math.ceil(math.log10(num_unique_chars))

# --- Step 4: Calculate the maximum number of characters (K) ---
# This is the total data available divided by the data needed for one character.
# We use integer division because we can only decode whole characters.
max_chars_k = total_data_size_d // digits_per_char

# --- Step 5: Print the final equation and result ---
print("Here is the calculation for K:")
print(f"K = (Total Symbols - (Number of Instructions * Size of one Instruction)) / Digits per Character")
print(f"K = ({total_symbols} - ({num_instructions} * {instruction_size_d})) / {digits_per_char}")
print(f"K = ({total_symbols} - {total_code_size_d}) / {digits_per_char}")
print(f"K = {total_data_size_d} / {digits_per_char}")
print(f"K = {max_chars_k}")

# Final Answer
print(f"<<<{max_chars_k}>>>")