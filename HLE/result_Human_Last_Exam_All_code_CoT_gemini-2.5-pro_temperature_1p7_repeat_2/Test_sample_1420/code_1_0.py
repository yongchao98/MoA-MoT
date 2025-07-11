# Define the parameters for the calculation
d = 3      # Spatial dimension
d_c = 4    # Upper critical dimension for the Ising universality class
n = 1      # Number of components of the order parameter (n=1 for the Ising model)

# Announce the task and the assumptions made
print("Calculating the specific heat scaling exponent α for d=3 using the epsilon expansion.")
print(f"We assume the system belongs to the Ising universality class (n = {n}).\n")

# --- Step 1: Calculate epsilon ---
# Epsilon (ϵ) is the small parameter for the expansion, defined as d_c - d
epsilon = d_c - d
print(f"Step 1: Calculate the expansion parameter ϵ.")
print(f"ϵ = d_c - d = {d_c} - {d} = {epsilon}\n")


# --- Step 2: Apply the formula for alpha ---
# The formula for α to first order in ϵ is α = (4 - n) / (2 * (n + 8)) * ϵ
print("Step 2: Use the first-order epsilon expansion formula for α.")
print("The formula is: α ≈ (4 - n) / (2 * (n + 8)) * ϵ\n")

# --- Step 3: Substitute values and calculate ---
# Substitute the numerical values into the formula and show the calculation
print("Step 3: Substitute the values and compute the result.")
print(f"Substituting n = {n} and ϵ = {epsilon}:")

numerator = 4 - n
denominator_term_1 = 2
denominator_term_2 = n + 8
full_denominator = denominator_term_1 * denominator_term_2

# Print the equation with all numbers filled in
print(f"α ≈ ({4} - {n}) / ({denominator_term_1} * ({n} + {8})) * {epsilon}")
print(f"α ≈ {numerator} / ({denominator_term_1} * {denominator_term_2}) * {epsilon}")
print(f"α ≈ {numerator} / {full_denominator} * {epsilon}")

# Final calculation
alpha = (numerator / full_denominator) * epsilon

# Print the final result in fractional and decimal form
print(f"\nFinal result as a fraction: α ≈ {numerator}/{full_denominator} = 1/6")
print(f"Final result as a decimal: α ≈ {alpha:.4f}")
