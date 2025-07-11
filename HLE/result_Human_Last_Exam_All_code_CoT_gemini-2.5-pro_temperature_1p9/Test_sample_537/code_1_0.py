# The pattern is the last digit of the square of consecutive integers n=1, 2, 3, 4, ...
# We need to find the 4th number in the sequence to decode the last word.
n = 4

# Calculate the square of n.
square = n * n

# Find the last digit of the square.
# The decoded number is this last digit.
decoded_number = square % 10

# Print the final equation step by step, as requested.
# This shows how the 4th number in the sequence is calculated.
print(f"The calculation for the 4th term is:")
print(n, "*", n, "=", square)
print(f"The decoded value of the last word is the last digit of the result.")
print(f"Decoded number: {decoded_number}")
