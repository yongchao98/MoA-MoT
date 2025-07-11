# The problem asks for a(3), where a(n) is the number of non-vanishing
# Feynman diagrams of order 2n.
# The formula is a(n) = (2n-1)!!

# Set the value for n
n = 3

# Calculate the starting number for the double factorial
start_num = 2 * n - 1

# Initialize the result for multiplication
result = 1
# Create a list to store the numbers for the equation string
factors = []

# Loop from the starting number down to 1 with a step of -2
for i in range(start_num, 0, -2):
    result *= i
    factors.append(str(i))

# Create the equation string by joining the factors with " * "
equation_str = " * ".join(factors)

# Print the final equation showing all the numbers and the result
print(f"{equation_str} = {result}")