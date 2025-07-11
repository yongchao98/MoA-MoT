import math

# Step 1: Calculate the total size of the program code.
num_instructions = 108
opcode_size = 1  # in decimal digits (D)
register_size = 1  # in decimal digits (D)
operand_size = 4  # in decimal digits (D)
instruction_size = opcode_size + register_size + operand_size
code_section_size = num_instructions * instruction_size

print(f"Each instruction has a size of {opcode_size}D (opcode) + {register_size}D (register) + {operand_size}D (operand) = {instruction_size}D.")
print(f"The program code consists of {num_instructions} instructions.")
print(f"Total size of the code section = {num_instructions} instructions * {instruction_size} D/instruction = {code_section_size} D.")
print("-" * 20)

# Step 2: Calculate the total size of the data.
total_symbols_in_ching = 9999
data_section_size = total_symbols_in_ching - code_section_size

print(f"The Ching contains a total of {total_symbols_in_ching} symbols (digits).")
print(f"Total size of the data section = {total_symbols_in_ching} D (total) - {code_section_size} D (code) = {data_section_size} D.")
print("-" * 20)

# Step 3: Determine the number of digits needed per character.
num_unique_chars = 108
# We need to find the smallest integer 'n' such that 10^n >= num_unique_chars.
# This is equivalent to ceil(log10(num_unique_chars)).
digits_per_char = math.ceil(math.log10(num_unique_chars))

print(f"There are {num_unique_chars} unique characters to represent.")
print("To represent all characters, we need to find the number of decimal digits 'n' where 10^n >= 108.")
print("10^2 = 100 (not enough)")
print("10^3 = 1000 (sufficient)")
print(f"Therefore, we need {digits_per_char} digits to represent one character.")
print("-" * 20)

# Step 4: Calculate the maximum number of characters (K).
# We use integer division since we can't decode a fraction of a character.
K = data_section_size // digits_per_char

print("The highest number of characters (K) that can be decoded is the total data size divided by the digits needed per character.")
print(f"K = {data_section_size} D / {digits_per_char} D/character = {K} characters.")

# The final answer K is the result of the calculation.
# print(f"\n<<<K = {K}>>>")
print(f"\n<<<{K}>>>")
