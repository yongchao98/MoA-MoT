# Define the transition rates from the problem description
lambda_01 = 0.019
lambda_10 = 0.65
lambda_12 = 0.4
lambda_21 = 0.392
lambda_23 = 0.008
lambda_31 = 0.008

# Explain the setup for finding the steady-state probabilities p_i
print("Step 1: Set up the steady-state equations by setting the derivatives to zero.")
print("The resulting linear equations are:")
print(f"  (1) 0 = -{lambda_01}*p0 + {lambda_10}*p1")
print(f"  (2) 0 = {lambda_12}*p1 - ({lambda_21} + {lambda_23})*p2")
print(f"  (3) 0 = {lambda_23}*p2 - {lambda_31}*p3")
print("  (4) p0 + p1 + p2 + p3 = 1 (Normalization condition)\n")

# Explain the process of solving the system
print("Step 2: Solve for p0, p2, and p3 in terms of p1.")
# Calculate the ratios
ratio_0_1 = lambda_10 / lambda_01
lambda_21_plus_23 = lambda_21 + lambda_23
ratio_2_1 = lambda_12 / lambda_21_plus_23
ratio_3_2 = lambda_23 / lambda_31
ratio_3_1 = ratio_3_2 * ratio_2_1

print(f"From (1), p0 = ({lambda_10}/{lambda_01})*p1 = {ratio_0_1:.4f}*p1")
print(f"From (2), p2 = ({lambda_12}/({lambda_21}+{lambda_23}))*p1 = ({lambda_12}/{lambda_21_plus_23})*p1 = {ratio_2_1:.4f}*p1")
print(f"From (3), p3 = ({lambda_23}/{lambda_31})*p2 = {ratio_3_2:.4f}*p2. Since p2 = {ratio_2_1:.4f}*p1, we get p3 = {ratio_3_1:.4f}*p1\n")

print("Step 3: Substitute these expressions into the normalization equation.")
print("p0 + p1 + p2 + p3 = 1")
print(f"({ratio_0_1:.4f}*p1) + p1 + ({ratio_2_1:.4f}*p1) + ({ratio_3_1:.4f}*p1) = 1")
denominator = ratio_0_1 + 1 + ratio_2_1 + ratio_3_1
print(f"p1 * ({ratio_0_1:.4f} + 1 + {ratio_2_1:.4f} + {ratio_3_1:.4f}) = 1")
print(f"p1 * {denominator:.4f} = 1 => p1 = 1 / {denominator:.4f}\n")


print("Step 4: Calculate the required sum P0(+inf) + P1(+inf) = p0 + p1.")
print("p0 + p1 = (ratio_0_1 * p1) + p1 = (ratio_0_1 + 1) * p1")
print("p0 + p1 = (ratio_0_1 + 1) / (ratio_0_1 + 1 + ratio_2_1 + ratio_3_1)\n")

print("Final equation with all numbers plugged in:")
numerator = ratio_0_1 + 1
# Print the final detailed equation as requested
equation = (
    f"p0 + p1 = (({lambda_10} / {lambda_01}) + 1) / "
    f"(({lambda_10} / {lambda_01}) + 1 + ({lambda_12} / ({lambda_21} + {lambda_23})) + "
    f"(({lambda_23} / {lambda_31}) * ({lambda_12} / ({lambda_21} + {lambda_23}))))"
)
print(equation)
# Simplify the equation with calculated ratios
simplified_equation = f"p0 + p1 = ({numerator:.4f}) / ({denominator:.4f})"
print(f"Which simplifies to:\n{simplified_equation}")


# Calculate the final result
result = numerator / denominator
print(f"\nThe final result for P0(+inf) + P1(+inf) is: {result}")