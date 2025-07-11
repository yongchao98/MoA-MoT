# Plan:
# 1. Define the parameters for the epsilon expansion: upper critical dimension (d_c), 
#    spatial dimension (d), and number of spin components (n for the Ising model).
# 2. Calculate epsilon (ϵ), which is the small parameter for the expansion, 
#    defined as ϵ = d_c - d.
# 3. Apply the first-order epsilon expansion formula for the specific heat exponent α for 
#    the n-component vector model: α = [(4 - n) / (2 * (n + 8))] * ϵ.
# 4. Substitute the numerical values into the formula to calculate α.
# 5. Print the full equation with all intermediate values to show the calculation process, 
#    followed by the final result.

# Parameters for the 3D Ising model
d_c = 4  # Upper critical dimension
d = 3    # Spatial dimension
n = 1    # Number of spin components (for the Ising model)

# Step 1: Calculate epsilon (ϵ)
epsilon = d_c - d

# Step 2: Calculate the numerator and denominator of the coefficient
numerator = 4 - n
denominator_full = 2 * (n + 8)

# Step 3: Calculate alpha (α) to first order in epsilon
alpha = (numerator / denominator_full) * epsilon

# Step 4: Print the calculation and result
print("The scaling exponent α for the specific heat is calculated using the first-order epsilon (ϵ) expansion.")
print("The general formula is: α ≈ (4 - n) / (2 * (n + 8)) * ϵ")
print(f"For a system with spatial dimension d = {d}, the upper critical dimension is d_c = {d_c}.")
print(f"Therefore, ϵ = d_c - d = {d_c} - {d} = {epsilon}.")
print(f"We use the Ising model, where the number of spin components is n = {n}.")
print("\nSubstituting these values into the formula:")
print(f"α ≈ ({4} - {n}) / (2 * ({n} + {8})) * {epsilon}")
print(f"α ≈ {numerator} / (2 * {n+8}) * {epsilon}")
print(f"α ≈ {numerator} / {denominator_full} * {epsilon}")
print(f"α ≈ {numerator/denominator_full} * {epsilon}")
print(f"α ≈ {alpha}")
print(f"As a decimal, α ≈ {alpha:.4f}")