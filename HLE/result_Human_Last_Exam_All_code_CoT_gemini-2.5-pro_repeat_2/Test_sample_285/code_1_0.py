import sympy

# We need to find the largest p for which the integral diverges.
# The condition for divergence is that the exponent in the radial integral is >= -1.
# The radial integral is of the form integral(R^k dR) from 1 to infinity.
# k = (dimension of subspace - 1) - p * (decay_rate)

# Dimension of the subspace of cubic coefficients (a6, a7, a8, a9) is 4.
dim_subspace = 4

# The decay rate of the integral I(a) in this subspace is ||a'||^(-2/3).
decay_rate = 2/3

# The exponent k in the radial integral is:
# k = (dim_subspace - 1) - p * decay_rate
# k = (4 - 1) - p * (2/3)
# k = 3 - 2*p/3

p = sympy.Symbol('p')
exponent = 3 - (2 * p / 3)

# The integral diverges if the exponent is >= -1.
divergence_condition = sympy.Ge(exponent, -1)

print("The condition for the integral to diverge is:")
print(f"{exponent} >= -1")

# Solve the inequality for p
solution = sympy.solve(divergence_condition, p)

print("\nSolving the inequality for p:")
# Step 1: Add 1 to both sides
step1 = sympy.Ge(exponent + 1, 0)
print(f"  {step1}")
# Step 2: Add 2*p/3 to both sides
step2 = sympy.Ge(4, 2 * p / 3)
print(f"  {step2}")
# Step 3: Multiply by 3
step3 = sympy.Ge(12, 2 * p)
print(f"  {step3}")
# Step 4: Divide by 2
step4 = sympy.Ge(6, p)
print(f"  {step4}")


print(f"\nThe function I(a) is not in L^p for p <= 6.")
# The largest p is the endpoint of this range.
largest_p = 6
print(f"The largest value of p such that I is not in L^p is {largest_p}.")

print("\nAs requested, the numbers in the final equation '3 - 2*p/3 >= -1' are:")
numbers = [3, 2, 3, -1]
for num in numbers:
    print(num)