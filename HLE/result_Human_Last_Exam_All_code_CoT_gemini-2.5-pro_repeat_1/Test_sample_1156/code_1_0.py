import math

# The normalised density of the invariant measure is rho(x) = 1 / (2 * ln(2) * (sqrt(x) + x))
# The code below will print this formula.

c1 = 1
c2 = 2
log_of_c2 = math.log(c2)

print("The normalised density of the invariant measure is given by the formula:")
# Output each number/component of the final equation
print(f"rho(x) = {c1} / ({c2} * ln({c2}) * (sqrt(x) + x))")
print(f"rho(x) = {c1} / ({c2} * {log_of_c2:.6f} * (sqrt(x) + x))")
final_denominator_const = c2 * log_of_c2
print(f"rho(x) = {c1} / ({final_denominator_const:.6f} * (sqrt(x) + x))")