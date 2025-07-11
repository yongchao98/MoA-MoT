import math

# Step 1: Identify the indices from the plots based on physical reasoning.
# n_0: base case simulation.
n_0 = 1
# n_1: 5x reduction in annulus pressure gradient (ΔP_a), slowing annulus flow. This matches plot 2.
n_1 = 2
# n_6: 5x reduction in viscosity ratio (μ = μ_tube/μ_annulus), increasing tube velocity. This matches plot 6.
n_6 = 6

# Step 2: Calculate the first term of the expression.
term1 = (n_1 / n_6) ** n_0

# Step 3: Determine the velocity ratio.
# For the given conditions, the ratio of maximum velocities simplifies.
# v_a_max / v_t_max = ln(3)
velocity_ratio = math.log(3)

# Step 4: Calculate the exponential term.
term2 = math.exp(velocity_ratio)

# Step 5: Calculate the final result.
final_result = term1 * term2

# Print the full equation with the calculated values
print(f"Based on the analysis of the plots, we have:")
print(f"n_0 = {n_0}")
print(f"n_1 = {n_1}")
print(f"n_6 = {n_6}")
print("\nThe first part of the expression is:")
print(f"({n_1}/{n_6})^{n_0} = {term1:.4f}")
print("\nFor the given physical conditions, the velocity ratio is:")
print(f"v_a,max / v_t,max = ln(3) ≈ {velocity_ratio:.4f}")
print("\nThe exponential term is:")
print(f"exp(v_a,max / v_t,max) = exp(ln(3)) = {term2:.4f}")
print("\nCombining these parts, the final calculation is:")
print(f"({n_1}/{n_6})^{n_0} * exp(v_a,max / v_t,max) = ({n_1}/{n_6})^{{{n_0}}} * {term2:.2f} = {term1:.4f} * {term2:.2f} = {final_result:.4f}")
print(f"\nThe final equation is ({n_1}/{n_6})^{n_0} * exp({velocity_ratio:.2f}) = {final_result:.2f}")
