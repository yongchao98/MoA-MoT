# Set the parameters for the vector space V = F_q^n
n = 3
q = 11

# The number of internal adjunctions is the order of the general linear group GL(n, q).
# This is calculated by finding the number of ways to choose n linearly independent vectors
# from an n-dimensional vector space over F_q.

print(f"The number of internal adjunctions from F_{q}^{n} to itself is the order of GL({n}, F_{q}).")
print("This is the number of invertible 3x3 matrices over F_11.")
print("The calculation involves multiplying the number of choices for each column of the matrix:")

# Calculate the number of choices for each column
q_n = q**n
choices_col1 = q_n - q**0
choices_col2 = q_n - q**1
choices_col3 = q_n - q**2

# Calculate the total number of adjunctions
total_adjunctions = choices_col1 * choices_col2 * choices_col3

# Print the calculation steps as requested
print(f"1. Choices for the first column (any non-zero vector): {q}^{n} - 1 = {choices_col1}")
print(f"2. Choices for the second column (not in the span of the first): {q}^{n} - {q} = {choices_col2}")
print(f"3. Choices for the third column (not in the span of the first two): {q}^{n} - {q**2} = {choices_col3}")

# Print the final equation with each number
print("\nThe final equation is the product of these choices:")
print(f"{choices_col1} * {choices_col2} * {choices_col3} = {total_adjunctions}")