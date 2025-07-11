import math

# The letters that can be arranged in any order are 'L', 'N', and 'S'.
# The letter 'W' must come last.
prefix_letters = ['L', 'N', 'S']

# The number of ways to arrange these letters is the factorial of the count of these letters.
n = len(prefix_letters)

# Calculate the factorial (number of arrangements).
result = math.factorial(n)

# Create the string for the equation, e.g., "3 * 2 * 1".
equation_parts = []
for i in range(n, 0, -1):
    equation_parts.append(str(i))
equation_str = " * ".join(equation_parts)

# Print the final equation.
print(f"The number of ways is the factorial of {n}:")
print(f"{equation_str} = {result}")
