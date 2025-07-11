import math

# The user can set the value of n here.
n = 15

# The maximum number of eigenvalues greater than 2 is given by the formula floor((n + 1) / 3).
# The final equation is: result = floor((n + 1) / 3)
# The numbers in this equation are n, 1, and 3.
val_n = n
val_1 = 1
val_3 = 3

# Calculate the numerator for the expression inside the floor function.
numerator = val_n + val_1

# Use integer division, which is equivalent to the floor function for positive numbers.
max_eigenvalues = numerator // val_3

print(f"To find the maximum number of eigenvalues greater than 2 for n = {val_n}, we use the formula:")
print("max_eigenvalues = floor((n + 1) / 3)")
print("\nBreaking down the calculation:")
print(f"n = {val_n}")
print(f"The numerator is n + 1 = {val_n} + {val_1} = {numerator}")
print(f"The denominator is {val_3}")
print(f"We need to compute the floor of {numerator} / {val_3}.")
print(f"The result is {max_eigenvalues}.")
