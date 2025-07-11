# Step 1: Define the constants from the problem description.
total_symbols_in_ching = 9999
num_instructions = 108
# Size of one instruction is 1D (opcode) + 1D (register) + 4D (operand).
instruction_size_in_digits = 6
# There are 108 possible characters to be decoded.
num_unique_characters = 108

# Step 2: Calculate the size of the program code.
program_code_size = num_instructions * instruction_size_in_digits

# Step 3: Calculate the size of the data available for decoding.
# This is the remainder of the Ching after the program code.
data_size = total_symbols_in_ching - program_code_size

# Step 4: Determine the number of digits required to represent one character.
# To represent one of 108 unique items (e.g., as numbers 0-107), we need 3 decimal digits.
digits_per_character = 3

# Step 5: Calculate K, the highest number of characters that can be decoded.
# This is the total data size divided by the data needed for one character.
K = data_size // digits_per_character

# Step 6: Print the result, showing the numbers in the final equation.
print(f"The total amount of data available is {data_size} decimal digits.")
print(f"Each character requires {digits_per_character} decimal digits to be represented.")
print("The highest number of characters (K) is calculated as:")
print(f"{data_size} / {digits_per_character} = {K}")