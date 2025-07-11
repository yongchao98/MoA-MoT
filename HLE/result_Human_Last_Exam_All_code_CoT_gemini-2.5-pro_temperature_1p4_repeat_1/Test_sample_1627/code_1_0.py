import fractions

# Step 1: Define the relationship between the radii.
# We are given that the outer radius 'b' is twice the inner radius 'a'.
# For calculation purposes, we can assume a = 1, which makes b = 2.
a = 1
b = 2 * a

# Step 2: The formula for the ratio of maximum tangential stress (sigma_t_max)
# to internal pressure (P_i) is (b^2 + a^2) / (b^2 - a^2).
# Let's calculate the numerator and denominator separately.
numerator_val = b**2 + a**2
denominator_val = b**2 - a**2

# Step 3: Calculate the final ratio.
ratio = numerator_val / denominator_val

# Step 4: Print the explanation and the step-by-step calculation.
print("The formula for the ratio of maximum tangential stress to internal pressure is:")
print("Ratio = (b^2 + a^2) / (b^2 - a^2)")
print("\nGiven that the outer radius 'b' is twice the inner radius 'a', we can set a=1 and b=2.")
print("Substituting these values into the formula:")
print(f"Ratio = ({b}^2 + {a}^2) / ({b}^2 - {a}^2)")
print(f"Ratio = ({int(b**2)} + {int(a**2)}) / ({int(b**2)} - {int(a**2)})")
print(f"Ratio = {int(numerator_val)} / {int(denominator_val)}")

# To show it as a simplified fraction.
ratio_fraction = fractions.Fraction(numerator_val, denominator_val)

print(f"\nThe calculated ratio is {ratio:.4f}, which is the fraction {ratio_fraction.numerator}/{ratio_fraction.denominator}.")
