# Define the parameters for the calculation
d = 3  # Spatial dimension
n = 1  # Number of components of the order parameter (for the Ising model)
d_c = 4 # Upper critical dimension for the n-vector model

# Calculate epsilon (ϵ)
epsilon = d_c - d

# Calculate the scaling exponent alpha (α) to first order in epsilon
# The formula is α = (4 - n) / (2 * (n + 8)) * ϵ
numerator = 4 - n
denominator = 2 * (n + 8)
alpha = (numerator / denominator) * epsilon

# Print the explanation and the final result
print("Calculating the specific heat exponent α using the epsilon expansion.")
print(f"Spatial dimension, d = {d}")
print(f"Upper critical dimension, d_c = {d_c}")
print(f"Number of order parameter components, n = {n} (Ising model)")
print("")
print("The expansion parameter epsilon (ϵ) is calculated as:")
print(f"ϵ = d_c - d = {d_c} - {d} = {epsilon}")
print("")
print("The first-order formula for α is: α = (4 - n) / (2 * (n + 8)) * ϵ")
print("Substituting the values:")
print(f"α = ({4} - {n}) / (2 * ({n} + {8})) * {epsilon}")
print(f"α = {numerator} / {denominator} * {epsilon}")
print(f"α = {numerator/denominator} * {epsilon}")
print("")
print(f"The final calculated value for α is: {alpha}")
print(f"As a fraction, α = 1/6")