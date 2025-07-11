# Parameters for the calculation
# d: spatial dimension
# d_c: upper critical dimension for the Ising universality class
# n: number of components of the order parameter (n=1 for the scalar Ising model)

d = 3
d_c = 4
n = 1

# Step 1: Calculate epsilon (ε)
epsilon = d_c - d

# Step 2: Calculate the numerator and denominator for the alpha formula
# The first-order epsilon expansion formula for alpha is:
# α = (4 - n) / (2 * (n + 8)) * ε
numerator = 4 - n
denominator = 2 * (n + 8)

# Step 3: Calculate alpha
alpha = (numerator / denominator) * epsilon

# Step 4: Print the explanation and the result, showing each number in the equation.
print("Calculating the specific heat scaling exponent α using the first-order epsilon expansion.")
print("-" * 75)
print(f"The given spatial dimension is d = {d}.")
print(f"The upper critical dimension for this system (Ising-like) is d_c = {d_c}.")
print(f"The expansion parameter is ε = d_c - d = {d_c} - {d} = {epsilon}.")
print(f"For the Ising model, the order parameter has n = {n} component.")
print("\nThe first-order formula for α is: α ≈ (4 - n) / (2 * (n + 8)) * ε")
print("\nSubstituting the values into the formula:")
print(f"α ≈ ({4} - {n}) / ({2} * ({n} + {8})) * {epsilon}")
print(f"α ≈ {numerator} / (2 * {n + 8}) * {epsilon}")
print(f"α ≈ {numerator} / {denominator} * {epsilon}")
print(f"α ≈ {alpha}")
print(f"\nThus, to first order in the epsilon expansion for d=3, the scaling exponent α is {alpha:.4f}.")
