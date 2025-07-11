import math

# Plan:
# 1. Define the geometric parameters of the system. The chord 'c' can be set to 1
#    as it will cancel out in the dimensionless factor K.
# 2. Calculate the dimensionless geometric interaction factor 'K' based on the formula derived
#    from the method of images and thin aerofoil theory:
#    K = (c/2) * [1/s - s / (s^2 + 4*h^2)]
# 3. Calculate the lift ratio L1/L2 using the formula:
#    L1/L2 = (1 + K) / (1 - K)
# 4. Print all intermediate values and the final result.

print("Step 1: Define geometric parameters")
# The chord length 'c' can be set to 1.0 for simplicity.
c = 1.0
# The separation 's' is 1/2 * c.
s = 0.5 * c
# The ride height 'h' is c/2.
h = 0.5 * c

print(f"Chord length, c = {c:.1f}")
print(f"Aerofoil separation, s = {s/c:.1f}c = {s:.1f}")
print(f"Ride height, h = {h/c:.1f}c = {h:.1f}")
print("-" * 40)

print("Step 2: Calculate the geometric interaction factor K")
# This factor represents the influence of the aerofoils on each other,
# including the ground effect.
term1_K = 1 / s
term2_K = s / (s**2 + 4 * h**2)
K = (c / 2) * (term1_K - term2_K)

print("The formula for K is: K = (c/2) * [1/s - s / (s^2 + 4*h^2)]")
print(f"Calculating the terms inside the bracket:")
print(f"  1/s = 1/{s:.1f} = {term1_K:.4f}")
print(f"  s / (s^2 + 4*h^2) = {s:.1f} / ({s:.1f}^2 + 4*{h:.1f}^2) = {term2_K:.4f}")
print(f"K = ({c:.1f}/2) * [{term1_K:.4f} - {term2_K:.4f}]")
print(f"K = {c/2:.1f} * [{term1_K - term2_K:.4f}]")
print(f"The calculated value of K is: {K:.4f}")
print("-" * 40)

print("Step 3: Calculate the final lift ratio L1/L2")
# The lift ratio L1/L2 is derived from solving the system of equations for the circulations.
numerator = 1 + K
denominator = 1 - K
lift_ratio = numerator / denominator

print("The formula for the lift ratio is: L1/L2 = (1 + K) / (1 - K)")
print(f"Plugging in the value of K:")
print(f"  Numerator = 1 + {K:.4f} = {numerator:.4f}")
print(f"  Denominator = 1 - {K:.4f} = {denominator:.4f}")
print(f"L1/L2 = {numerator:.4f} / {denominator:.4f}")

# Check if the result is a whole number for cleaner printing
if lift_ratio == int(lift_ratio):
    lift_ratio = int(lift_ratio)

print("-" * 40)
print(f"The final calculated lift ratio L1/L2 is: {lift_ratio}")