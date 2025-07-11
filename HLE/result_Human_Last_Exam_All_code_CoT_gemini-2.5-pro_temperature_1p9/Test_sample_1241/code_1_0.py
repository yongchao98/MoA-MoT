import numpy as np

# This script calculates the steady-state probability sum P0 + P1.

# Plan:
# 1. Use the steady-state conditions (P_i' = 0) from the (corrected)
#    Kolmogorov-Chapman equations to establish relationships between the
#    probabilities P0, P1, P2, and P3.
# 2. Use the given numerical values for lambda to simplify these relationships.
# 3. Apply the normalization condition (P0 + P1 + P2 + P3 = 1) to find the
#    explicit values for the probabilities.
# 4. Compute the final sum P0 + P1 and display the full calculation.

# Given lambda values
lambda_01 = 0.019
lambda_10 = 0.65
lambda_12 = 0.4
lambda_21 = 0.392
lambda_23 = 0.008
lambda_31 = 0.008

# From the steady-state equations, we derive relations between probabilities.
# Equation 3: lambda_12*P1 = (lambda_21 + lambda_23)*P2
# Equation 4: lambda_23*P2 = lambda_31*P3
# Equation 1: lambda_01*P0 = lambda_10*P1

print("We are solving the steady-state equations for P0, P1, P2, and P3.")
print("The values of the transition rates are:")
print(f"lambda_01 = {lambda_01}, lambda_10 = {lambda_10}, lambda_12 = {lambda_12}, lambda_21 = {lambda_21}, lambda_23 = {lambda_23}, lambda_31 = {lambda_31}\n")

print("From the equations, we find relationships between the probabilities by substituting the given values:")
# Check P1, P2, P3 relationship
p1_p2_coeff = (lambda_21 + lambda_23) / lambda_12
p2_p3_coeff = lambda_31 / lambda_23
print(f"  {lambda_12}*P1 = ({lambda_21} + {lambda_23})*P2  => {lambda_12}*P1 = {lambda_21+lambda_23:.3f}*P2. Since {lambda_12} = {lambda_21+lambda_23:.3f}, we have P1 = P2.")
print(f"  {lambda_23}*P2 = {lambda_31}*P3  => {lambda_23}*P2 = {lambda_31}*P3. Since {lambda_23} = {lambda_31}, we have P2 = P3.")
print("Therefore, P1 = P2 = P3.\n")

# P0 relation
p0_p1_coeff = lambda_10 / lambda_01
print(f"Also, from {lambda_01}*P0 = {lambda_10}*P1, we get P0 = ({lambda_10}/{lambda_01})*P1 = {p0_p1_coeff:.5f} * P1\n")


# Normalization
# P0 + P1 + P2 + P3 = 1  => (p0_p1_coeff)*P1 + P1 + P1 + P1 = 1
# (p0_p1_coeff + 3)*P1 = 1
total_coeff_p1 = p0_p1_coeff + 3
p1 = 1 / total_coeff_p1
p0 = p0_p1_coeff * p1

print("Using the normalization condition P0 + P1 + P2 + P3 = 1:")
print(f"  ({p0_p1_coeff:.5f})*P1 + P1 + P1 + P1 = 1")
print(f"  ({total_coeff_p1:.5f})*P1 = 1 => P1 = {p1:.5f}")
print(f"This gives us the individual probabilities:")
print(f"  P0 = {p0:.5f}, P1 = {p1:.5f}, P2 = {p1:.5f}, P3 = {p1:.5f}\n")


# Final calculation
# The sum P0 + P1 can be expressed as a function of the lambdas:
# P0 + P1 = (lambda_10/lambda_01 + 1) * P1
# P0 + P1 = (lambda_10/lambda_01 + 1) / (lambda_10/lambda_01 + 3)
# P0 + P1 = (lambda_10 + lambda_01) / (lambda_10 + 3*lambda_01)
numerator = lambda_10 + lambda_01
denominator = lambda_10 + 3 * lambda_01
final_result = numerator / denominator

print("The quantity to find is P0 + P1.")
print("The final calculation is expressed as an equation with the given parameters:")
print(f"P0 + P1 = (lambda_10 + lambda_01) / (lambda_10 + 3 * lambda_01)")
print(f"Plugging in the numbers:")
print(f"P0 + P1 = ({lambda_10} + {lambda_01}) / ({lambda_10} + 3 * {lambda_01})")
print(f"P0 + P1 = {numerator} / {denominator}")
print(f"P0 + P1 = {final_result}")