# The provided Javascript is obfuscated and buggy, but contains a clear clue.
# It constructs the number 9 by adding 'true' (coerced to 1) nine times.
# This script builds and prints the equation that demonstrates this calculation,
# as requested by the prompt.

# The list of numbers in the equation
numbers = [1] * 9

# The final result
result = sum(numbers)

# Create the equation string like "1 + 1 + ... + 1"
equation_string = " + ".join(map(str, numbers))

# Print the final equation, showing each number, as requested.
print(f"{equation_string} = {result}")