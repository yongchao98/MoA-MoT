import math

# The least upper bound is achieved at n=1.
n = 1

# Calculate the four terms for n=1.
# E_H = (n+2)/(n+1)
E_H = (n + 2) / (n + 1)

# S_H = n+1
S_H = float(n + 1)

# S_P = 1, as P is orthogonal
S_P = 1.0

# E_P <= 2/(n+1). This bound is achieved for n=1, so E_P = 2/(n+1).
E_P = 2 / (n + 1)

# Calculate the final product
product = E_P * E_H * S_P * S_H

# Print the equation with the calculated values
print(f"The maximum value is achieved for n = {n}.")
print("The product is given by the equation: E_P * E_H * S_P * S_H")
print(f"The values are: {E_P} * {E_H} * {S_P} * {S_H} = {product}")
print("The least upper bound of the product is the value calculated for n=1.")

# The final answer in the required format
# print(f"<<<{product}>>>")