import sys

# Set a higher precision for floating-point output
# This is optional but can be useful for verification
# sys.float_info.dig = 15

# Step 1: Define the given transition rates (lambda values)
lambda_01 = 0.019
lambda_10 = 0.65
lambda_12 = 0.4
lambda_21 = 0.392
lambda_23 = 0.008
lambda_31 = 0.008

# Step 2: Express P0, P2, and P3 in terms of P1
# From P0'(t) = 0: -lambda_01*P0 + lambda_10*P1 = 0  => P0 = (lambda_10 / lambda_01) * P1
# Let's call the coefficient c0
c0 = lambda_10 / lambda_01

# From P2'(t) = 0: lambda_12*P1 - (lambda_21 + lambda_23)*P2 = 0 => P2 = (lambda_12 / (lambda_21 + lambda_23)) * P1
# Let's call the coefficient c2
c2 = lambda_12 / (lambda_21 + lambda_23)

# From P3'(t) = 0: lambda_23*P2 - lambda_31*P3 = 0 => P3 = (lambda_23 / lambda_31) * P2
# Substitute P2: P3 = (lambda_23 / lambda_31) * c2 * P1
# Let's call the coefficient c3
c3 = (lambda_23 / lambda_31) * c2

# Step 3: Use the normalization condition P0 + P1 + P2 + P3 = 1
# c0*P1 + P1 + c2*P1 + c3*P1 = 1
# P1 * (c0 + 1 + c2 + c3) = 1
# P1 = 1 / (c0 + 1 + c2 + c3)

# We need to find P0 + P1
# P0 + P1 = c0*P1 + P1 = (c0 + 1) * P1
# P0 + P1 = (c0 + 1) / (c0 + 1 + c2 + c3)

# Step 4: Calculate the final result
numerator = c0 + 1
denominator = c0 + 1 + c2 + c3
result = numerator / denominator

# Step 5: Print the final equation with the numbers and the result
print("The relationship between the steady-state probabilities can be expressed in terms of P1:")
print(f"P0 = ({lambda_10}/{lambda_01}) * P1 = {c0:.4f} * P1")
print(f"P2 = ({lambda_12}/({lambda_21}+{lambda_23})) * P1 = {c2:.4f} * P1")
print(f"P3 = ({lambda_23}/{lambda_31}) * P2 = {c3:.4f} * P1\n")

print("The quantity P0(inf) + P1(inf) is calculated as (c0 + 1) / (c0 + 1 + c2 + c3).")
print("Plugging in the numbers for the coefficients:")
print(f"P0(inf) + P1(inf) = ({c0:.15f} + 1) / ({c0:.15f} + 1 + {c2:.1f} + {c3:.1f})")
print(f"P0(inf) + P1(inf) = {numerator:.15f} / {denominator:.15f}")
print(f"P0(inf) + P1(inf) = {result:.15f}")