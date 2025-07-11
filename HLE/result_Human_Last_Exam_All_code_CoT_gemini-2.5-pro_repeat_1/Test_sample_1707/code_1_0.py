# Define the initial values as strings for clear symbolic representation
alpha_str = "10**10000"
x0_str = "10**(-5000000)"

# The formula for T is derived from the solvability condition of the first-order perturbation problem.
# The derivation leads to the equation: T = alpha * (1 - x0) / (2 * x0).

print("The derived formula for T is:")
print("T = alpha * (1 - x0) / (2 * x0)")

print("\nSubstituting the given values:")
print(f"T = ({alpha_str}) * (1 - {x0_str}) / (2 * {x0_str})")

print("\nSimplifying the expression by splitting the fraction:")
# T = (alpha / (2*x0)) * (1 - x0)
# First, calculate the term alpha / (2 * x0)
# (10**10000) / (2 * 10**-5000000) = (1/2) * 10**(10000 - (-5000000)) = 0.5 * 10**5010000
coeff_part1 = 0.5
exp_part1 = 10000 - (-5000000)
print(f"The term (alpha / (2 * x0)) simplifies to: {coeff_part1} * 10**{exp_part1}")

print("\nSo the equation for T becomes:")
print(f"T = ({coeff_part1} * 10**{exp_part1}) * (1 - {x0_str})")

print("\nFinally, we distribute the term to get the final equation:")
# T = (0.5 * 10**5010000) * 1 - (0.5 * 10**5010000) * (10**-5000000)
# T = 0.5 * 10**5010000 - 0.5 * 10**(5010000 - 5000000)
# T = 0.5 * 10**5010000 - 0.5 * 10**10000
coeff_part2 = coeff_part1
exp_part2 = exp_part1 - 5000000

print(f"T = {coeff_part1} * 10**{exp_part1} - {coeff_part1} * 10**({exp_part1} - 5000000)")
print("\nWhich results in the final symbolic form:")
print(f"T = {coeff_part1} * 10**{exp_part1} - {coeff_part2} * 10**{exp_part2}")