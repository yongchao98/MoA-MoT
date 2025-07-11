import math

# --- Step 1: Calculate the total size of the program code ---
ching_total_symbols = 9999
num_instructions = 108

# Each instruction has format [opcode][register][operand]
# opcode is 1 digit, register is 1 digit, operand is 4 digits
instruction_size_d = 1 + 1 + 4
code_size_d = num_instructions * instruction_size_d

print(f"The Ching contains {ching_total_symbols} symbols (decimal digits).")
print(f"The program has {num_instructions} instructions.")
print(f"The size of one instruction is {instruction_size_d} decimal digits (D).")
print(f"Total code size = {num_instructions} instructions * {instruction_size_d} D/instruction = {code_size_d} D.")
print("-" * 20)

# --- Step 2: Calculate the available space for data ---
data_size_d = ching_total_symbols - code_size_d

print(f"Total data size = {ching_total_symbols} D (Ching size) - {code_size_d} D (code size) = {data_size_d} D.")
print("-" * 20)

# --- Step 3: Determine the minimum digits required per character ---
num_unique_chars = 108

# We need to find the smallest integer 'd' such that 10^d >= num_unique_chars
# This is equivalent to ceil(log10(num_unique_chars))
digits_per_char = math.ceil(math.log10(num_unique_chars))

print(f"There are {num_unique_chars} unique characters to encode.")
print("To represent each character, we need a number of decimal digits 'd' where 10^d >= 108.")
print(f"10^2 = 100 (not enough)")
print(f"10^3 = 1000 (enough)")
print(f"So, the minimum digits needed per character is {digits_per_char} D.")
print("-" * 20)

# --- Step 4: Calculate the highest number of characters (K) ---
# We use integer division since we can't decode a fraction of a character.
max_chars_k = data_size_d // digits_per_char

print("The highest number of characters (K) is the total data size divided by the data size per character.")
print(f"K = {data_size_d} D / {digits_per_char} D/character = {max_chars_k} characters.")
print("\nTherefore, the final answer for K is:")
print(max_chars_k)
