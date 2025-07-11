import math

# Step 1: Calculate the total size of the program code.
num_instructions = 108
# Each instruction has a fixed format: [opcode][register][operand]
# opcode = 1 Decimal Digit (D)
# register = 1 Decimal Digit (D)
# operand = 4 Decimal Digits (D)
instruction_size_d = 1 + 1 + 4
code_size_d = num_instructions * instruction_size_d

# Step 2: Calculate the available data size.
total_ching_size_d = 9999
data_size_d = total_ching_size_d - code_size_d

# Step 3: Determine the space needed for one character.
num_possible_characters = 108
# We need to find the smallest integer 'n' such that 10^n >= 108.
# This is equivalent to ceiling(log10(108)).
digits_per_character = math.ceil(math.log10(num_possible_characters))

# Step 4: Calculate the maximum number of characters (K).
# This is the total data size divided by the digits needed for one character.
# We use integer division since we can't have a fraction of a character.
K = data_size_d // digits_per_character

# Print the breakdown of the final calculation.
print(f"The program code occupies {num_instructions} * {instruction_size_d} = {code_size_d} digits.")
print(f"The available data space is {total_ching_size_d} - {code_size_d} = {data_size_d} digits.")
print(f"Each character requires {digits_per_character} digits to be uniquely represented (ceil(log10({num_possible_characters}))).")
print("\nThe final calculation for K, the highest number of characters, is:")
print(f"K = Data Space / Digits per Character")
print(f"K = {data_size_d} / {digits_per_character}")
print(f"K = {K}")