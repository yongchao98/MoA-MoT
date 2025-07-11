# Step 1: Define the constants based on the problem description.
total_digits_in_ching = 9999
num_instructions = 108
num_unique_characters = 108

# Size of each part of an instruction in decimal digits (D)
opcode_size = 1
register_size = 1
operand_size = 4

# Step 2: Calculate the size of a single instruction.
# Instruction format: [opcode][register][operand]
instruction_size = opcode_size + register_size + operand_size

# Step 3: Calculate the total size of the code section.
# The program has 108 instructions.
code_section_size = num_instructions * instruction_size

# Step 4: Calculate the size of the data section.
# The remaining space in the Ching is used for data.
data_section_size = total_digits_in_ching - code_section_size

# Step 5: Determine the number of digits required to encode one character.
# To represent 108 unique characters, we need to find the smallest power of 10
# that is greater than or equal to 108.
# 10^1 = 10 (not enough)
# 10^2 = 100 (not enough)
# 10^3 = 1000 (enough)
# So, we need 3 digits per character.
digits_per_character = 3

# Step 6: Calculate K, the highest number of characters that can be decoded.
# This is the total data size divided by the size of one character's encoding.
K = data_section_size // digits_per_character

# Final Step: Print the equation with all the numbers and the final answer.
# This shows the full calculation from start to finish.
print(f"Total data space = {total_digits_in_ching} (total digits) - ({num_instructions} (instructions) * {instruction_size} (digits/instruction)) = {data_section_size} digits")
print(f"Max characters (K) = {data_section_size} (data digits) / {digits_per_character} (digits/character) = {K}")
print("\nFinal equation:")
print(f"K = ({total_digits_in_ching} - ({num_instructions} * {instruction_size})) / {digits_per_character} = {K}")

<<<3117>>>