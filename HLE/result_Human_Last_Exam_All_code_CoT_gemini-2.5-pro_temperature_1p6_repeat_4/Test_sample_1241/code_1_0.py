# Define the given transition rates
lambda_01 = 0.019
lambda_10 = 0.65
lambda_12 = 0.4
lambda_21 = 0.392
lambda_23 = 0.008
lambda_31 = 0.008

# At steady state, the derivatives are zero, leading to a system of linear equations.
# As derived in the explanation, the solution for P0(+inf) + P1(+inf),
# which we denote as p0 + p1, is given by the formula:
# p0 + p1 = ((lambda_10 / lambda_01) + 1) / ((lambda_10 / lambda_01) + 3)

# Perform the calculation
ratio_10_01 = lambda_10 / lambda_01
numerator = ratio_10_01 + 1
denominator = ratio_10_01 + 3
result = numerator / denominator

# Print the final equation with the numerical values
print("The final equation for P0(+inf) + P1(+inf) is:")
print(f"P0(+inf) + P1(+inf) = (({lambda_10} / {lambda_01}) + 1) / (({lambda_10} / {lambda_01}) + 3)")
print("\nSubstituting the values:")
print(f"P0(+inf) + P1(+inf) = ({ratio_10_01:.5f} + 1) / ({ratio_10_01:.5f} + 3)")
print(f"P0(+inf) + P1(+inf) = {numerator:.5f} / {denominator:.5f}")
print(f"\nThe final result is: {result:.7f}")