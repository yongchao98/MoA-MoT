import fractions

# Step 1: Define the given ratio of outer radius to inner radius.
k = 2.0

# Step 2: The formula for the ratio of maximum tangential stress (sigma_t_max)
# to internal pressure (P_i) is (k^2 + 1) / (k^2 - 1), where k = r_o / r_i.
print("The formula for the ratio is: σ_t_max / P_i = (k^2 + 1) / (k^2 - 1)")
print(f"Given the outer radius is twice the inner radius, the ratio k = {int(k)}.")
print("\nSubstitute k = 2 into the formula:")

# Step 3: Calculate the numerator and denominator.
numerator = k**2 + 1
denominator = k**2 - 1

# Step 4: Display the calculation steps with the numbers.
print(f"σ_t_max / P_i = (({int(k)})^2 + 1) / (({int(k)})^2 - 1)")
print(f"σ_t_max / P_i = ({int(k**2)} + 1) / ({int(k**2)} - 1)")
print(f"σ_t_max / P_i = {int(numerator)} / {int(denominator)}")

# Step 5: Calculate the final ratio.
ratio_decimal = numerator / denominator
ratio_fraction = fractions.Fraction(numerator, denominator)

print(f"\nThe final ratio is {ratio_fraction}, which is approximately {ratio_decimal:.4f}.")
