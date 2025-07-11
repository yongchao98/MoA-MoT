# The Betti numbers for the Möbius band (which has the homotopy type of a circle).
# beta_0 is the number of connected components.
beta_0 = 1
# beta_1 is the number of 1-dimensional holes.
beta_1 = 1

# The sum of all Betti numbers (higher Betti numbers are zero).
sum_of_betti_numbers = beta_0 + beta_1

# The number of non-collapsible rooted forests is the square of the sum of the Betti numbers.
result = sum_of_betti_numbers**2

# Print the final calculation, showing each number in the equation.
print(f"The calculation is based on the Betti numbers of the Möbius band (β_0, β_1, ...).")
print(f"β_0 = {beta_0}")
print(f"β_1 = {beta_1}")
print(f"The number of non-collapsible forests is (β_0 + β_1)^2.")
print(f"So, the final equation is: ({beta_0} + {beta_1})^2 = {result}")
