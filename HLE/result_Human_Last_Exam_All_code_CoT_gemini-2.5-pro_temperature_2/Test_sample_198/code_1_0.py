# Step 1: Define the constants based on the problem description.
total_digits_in_ching = 9999
num_instructions = 108
num_unique_characters = 108

# Size of one instruction in decimal digits (D)
# Format: [opcode][register][operand]
opcode_size = 1
register_size = 1
operand_size = 4
instruction_size = opcode_size + register_size + operand_size

# Step 2: Calculate the total size of the program code.
program_code_size = num_instructions * instruction_size

# Step 3: Calculate the total size of the data segment.
data_size = total_digits_in_ching - program_code_size

# Step 4: Determine the number of digits needed to represent one character.
# We need to find the smallest integer 'd' such that 10^d >= num_unique_characters.
# 10^1 = 10 < 108
# 10^2 = 100 < 108
# 10^3 = 1000 > 108
# So, 3 digits are required per character.
digits_per_character = 3

# Step 5: Calculate K, the highest number of characters that can be decoded.
K = data_size // digits_per_character

# Final output: Print the equation and the result.
print(f"The total data size is {data_size} digits.")
print(f"Each character requires {digits_per_character} digits to be encoded.")
print("The final equation for K is: data_size / digits_per_character")
print(f"K = {data_size} / {digits_per_character} = {K}")
print("\nThe highest number of characters that can be decoded from the Ching is:")
print(K)
