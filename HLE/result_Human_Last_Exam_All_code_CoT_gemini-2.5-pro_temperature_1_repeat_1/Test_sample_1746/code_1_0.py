# Statement G states that for positive integers n, m, the identity
# (n // m) * m + n % m == n always holds true. This is a fundamental
# property of Python's integer division and modulo operators.

# Let's demonstrate this with an example.
n = 127
m = 12

# Calculate the quotient and remainder
quotient = n // m
remainder = n % m

# Reconstruct n using the identity
reconstructed_n = quotient * m + remainder

print(f"Demonstrating statement G for n = {n} and m = {m}:")
print(f"The identity is: (n // m) * m + n % m = n")
print(f"Calculating the components:")
print(f"Quotient (n // m): {quotient}")
print(f"Remainder (n % m): {remainder}")
print("\nSubstituting these values into the identity's final equation:")
# The following line prints each number in the final equation, as requested.
print(f"{quotient} * {m} + {remainder} = {n}")

# Verify the calculation
print(f"\nFinal check: {quotient * m + remainder} = {n}, which is {reconstructed_n == n}.")
