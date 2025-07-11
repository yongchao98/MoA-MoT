from fractions import Fraction

# Parameters for the calculation
d = 3      # Spatial dimension
d_c = 4    # Upper critical dimension for the O(N) model
N = 1      # Number of spin components (assuming the Ising model)

# --- Explanation ---
print("This script calculates the specific heat scaling exponent α for d=3 using the epsilon expansion from renormalization group theory.")
print("The general formula for α to the first order in epsilon (ϵ = d_c - d) for the O(N) model is:")
print("α = (4 - N) / (2 * (N + 8)) * ϵ\n")

# --- Step-by-step Calculation ---
print("Step 1: Define the system parameters.")
print(f"  - Spatial dimension, d = {d}")
print(f"  - Upper critical dimension, d_c = {d_c}")
print(f"  - Number of spin components, N = {N} (assuming the Ising model)\n")

print("Step 2: Calculate the expansion parameter ϵ.")
epsilon = d_c - d
print(f"  ϵ = d_c - d")
print(f"  ϵ = {d_c} - {d} = {epsilon}\n")

print("Step 3: Substitute the parameters into the formula for α.")
# Retrieve each number for the equation output
numerator_term_1 = 4
numerator_term_2 = N
denominator_term_1 = 2
denominator_term_2 = N
denominator_term_3 = 8
final_epsilon = epsilon

# Print the equation with all numbers
print(f"  α = ({numerator_term_1} - {numerator_term_2}) / ({denominator_term_1} * ({denominator_term_2} + {denominator_term_3})) * {final_epsilon}")

# Calculate intermediate values for clarity
numerator_val = numerator_term_1 - numerator_term_2
denominator_intermediate = denominator_term_2 + denominator_term_3
denominator_val = denominator_term_1 * denominator_intermediate
print(f"  α = {numerator_val} / ({denominator_term_1} * {denominator_intermediate}) * {final_epsilon}")
print(f"  α = {numerator_val} / {denominator_val} * {final_epsilon}\n")

print("Step 4: Compute the final value of α.")
alpha_fraction = Fraction(numerator_val * epsilon, denominator_val)
alpha_decimal = float(alpha_fraction)

print(f"  The exact value of α as a fraction is: {alpha_fraction.numerator}/{alpha_fraction.denominator}")
print(f"  The decimal value is α ≈ {alpha_decimal:.4f}")