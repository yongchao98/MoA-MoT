# Define the coefficients and powers
a = 7
b = -72
n = 8

# Calculate the antiderivative
antiderivative_constant = a
antiderivative_x_term = b / (n + 1)

# Print the result
print(f"{antiderivative_constant}*X + {antiderivative_x_term}*X**{n+1} + C")