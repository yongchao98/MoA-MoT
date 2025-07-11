import numpy as np

# Given transition rates
lambda_01 = 0.019
lambda_10 = 0.65
lambda_12 = 0.4
lambda_21 = 0.392
lambda_23 = 0.008
lambda_31 = 0.008

# At steady state, the derivatives are zero. We solve the resulting system of linear equations.
# P_0'(t) = 0  =>  lambda_01 * p0 = lambda_10 * p1
# P_2'(t) = 0  => (lambda_21 + lambda_23) * p2 = lambda_12 * p1
# P_3'(t) = 0  =>  lambda_31 * p3 = lambda_23 * p2
# And the normalization condition: p0 + p1 + p2 + p3 = 1

# Step 1: Express p0, p2, and p3 in terms of p1.
# From P_0'(t) = 0:
# p0 = (lambda_10 / lambda_01) * p1
coeff_p0 = lambda_10 / lambda_01

# From P_2'(t) = 0:
# p2 = (lambda_12 / (lambda_21 + lambda_23)) * p1
coeff_p2 = lambda_12 / (lambda_21 + lambda_23)

# From P_3'(t) = 0:
# p3 = (lambda_23 / lambda_31) * p2
# Substitute p2: p3 = (lambda_23 / lambda_31) * coeff_p2 * p1
coeff_p3 = (lambda_23 / lambda_31) * coeff_p2

# Step 2: Substitute into the normalization equation to solve for p1.
# p0 + p1 + p2 + p3 = 1
# (coeff_p0 * p1) + p1 + (coeff_p2 * p1) + (coeff_p3 * p1) = 1
# p1 * (coeff_p0 + 1 + coeff_p2 + coeff_p3) = 1
p1 = 1 / (coeff_p0 + 1 + coeff_p2 + coeff_p3)

# Step 3: Calculate p0 using the value of p1.
p0 = coeff_p0 * p1

# Step 4: Calculate the required sum p0 + p1.
result = p0 + p1

# Print the results as requested.
print(f"The steady-state probabilities are related as follows:")
print(f"p0 = ({lambda_10}/{lambda_01}) * p1 = {coeff_p0:.4f} * p1")
print(f"p2 = ({lambda_12}/({lambda_21}+{lambda_23})) * p1 = {coeff_p2:.4f} * p1")
print(f"p3 = (({lambda_23}/{lambda_31}) * {coeff_p2:.4f}) * p1 = {coeff_p3:.4f} * p1")
print("\nSolving for the probabilities:")
print(f"p0 = {p0:.4f}")
print(f"p1 = {p1:.4f}")
# For completeness, let's also calculate and print p2 and p3
p2 = coeff_p2 * p1
p3 = coeff_p3 * p1
print(f"p2 = {p2:.4f}")
print(f"p3 = {p3:.4f}")
print(f"Check sum: p0+p1+p2+p3 = {p0+p1+p2+p3:.4f}")


print("\nThe final equation is:")
print(f"P0(inf) + P1(inf) = {p0:.4f} + {p1:.4f} = {result:.4f}")
print(f"\nFinal result as a fraction: {int((p0+p1)*707)}/707")

# Store the precise final answer for submission
final_answer = (lambda_10 + lambda_01) / (lambda_10 + lambda_01 + lambda_12 + (lambda_12*lambda_23/lambda_31))

print(f"\nExact value: {final_answer}")