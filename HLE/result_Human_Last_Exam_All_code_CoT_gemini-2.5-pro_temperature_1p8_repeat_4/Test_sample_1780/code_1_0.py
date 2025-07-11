import math

# The problem reduces to finding the number of permutations of the 3 letters
# that can precede the terminal letter 'W'.
num_letters_to_arrange = 3

# The number of ways is n!, where n is the number of letters.
result = math.factorial(num_letters_to_arrange)

# Create the equation string for the factorial calculation.
# For n=3, this will be "3 * 2 * 1".
equation_parts = []
for i in range(num_letters_to_arrange, 0, -1):
    equation_parts.append(str(i))
equation_str = " * ".join(equation_parts)

# Print the final result along with the equation used for the calculation.
print(f"The number of arrangements is the number of permutations of 3 letters (L, N, S), which is 3!.")
print(f"The final calculation is: {equation_str} = {result}")
