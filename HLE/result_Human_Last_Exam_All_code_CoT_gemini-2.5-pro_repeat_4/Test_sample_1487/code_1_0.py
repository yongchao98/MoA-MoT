# The analytical derivation simplifies the problem to a simple arithmetic calculation.
# The expression simplifies to: 1 + 10^15

# Define the two terms of the final sum.
first_term = 1
second_term = 10**15

# Calculate the final result.
result = first_term + second_term

# As requested, print each number in the final equation.
# We use int() to avoid scientific notation for the large number.
print(f"{first_term} + {int(second_term)} = {int(result)}")