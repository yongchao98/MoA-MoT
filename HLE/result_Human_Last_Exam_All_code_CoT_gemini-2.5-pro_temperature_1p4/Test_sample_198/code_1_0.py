import math

# Step 1: Calculate the memory space occupied by the program code.
# The program has 108 instructions.
num_instructions = 108
# Each instruction has the format [opcode][register][operand], which is
# 1 digit + 1 digit + 4 digits = 6 digits long.
instruction_size_d = 6
code_space_d = num_instructions * instruction_size_d

# Step 2: Calculate the remaining memory space available for data.
# The Ching contains a total of 9999 symbols (digits).
total_ching_symbols_d = 9999
# The data is what remains after the code.
data_space_d = total_ching_symbols_d - code_space_d

# Step 3: Determine the minimum number of decimal digits required to represent a single character.
# The decoded message uses characters from a set of 108 unique characters.
num_unique_chars = 108
# To represent one of 108 unique items using decimal digits, we need to find the smallest
# integer 'n' such that 10^n is greater than or equal to 108.
# This can be calculated as the ceiling of log base 10 of 108.
digits_per_char = math.ceil(math.log10(num_unique_chars))

# Step 4: Calculate K, the highest number of characters that can be encoded.
# The highest number of characters is the total number of data digits divided by the
# number of digits needed to represent one character. We take the floor of the result
# because a partial encoding cannot produce a full character.
K = data_space_d // digits_per_char

# Final output explaining the calculation.
print(f"The highest number of characters (K) can be calculated as follows:")
print(f"1. Calculate total digits for code: {num_instructions} instructions * {instruction_size_d} digits/instruction = {code_space_d} digits.")
print(f"2. Calculate total digits for data: {total_ching_symbols_d} total digits - {code_space_d} code digits = {data_space_d} data digits.")
print(f"3. Digits to represent 1 of {num_unique_chars} characters: ceil(log10({num_unique_chars})) = {digits_per_char} digits.")
print(f"4. Calculate K by dividing data digits by digits per character:")
print(f"K = floor( (Total Digits - Code Digits) / Digits per Character )")
print(f"K = floor( ({total_ching_symbols_d} - {code_space_d}) / {digits_per_char} )")
print(f"K = floor( {data_space_d} / {digits_per_char} )")
print(f"K = {K}")
<<<3117>>>