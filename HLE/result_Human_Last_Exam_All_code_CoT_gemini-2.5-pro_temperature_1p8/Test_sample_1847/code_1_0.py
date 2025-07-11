import math

# This script calculates and prints the analytical solution for the specified shock wave density profile.

# Header for the solution
print("The analytical solution for the normalized density profile `rho(z) / rho_0` is given by the following equation:")
print("")  # For a new line

# Define the components of the solution to satisfy the prompt's request
# to output each number in the final equation.
eq_lhs = "rho(z) / rho_0"
num = 8.0
den_c1 = 8.0
exp1_coeff_num = 4.0
exp1_coeff_den = 3.0
exp2_coeff_num = 8.0
exp2_coeff_den = 3.0
sqrt_coeff = 8.0

exp1_coeff_str = f"({int(exp1_coeff_num)}/{int(exp1_coeff_den)})"
exp2_coeff_str = f"({int(exp2_coeff_num)}/{int(exp2_coeff_den)})"


# Print the final equation using the defined components
print(f"{eq_lhs} = {int(num)} / ({int(den_c1)} + exp({exp1_coeff_str}*z) - sqrt(exp({exp2_coeff_str}*z) + {int(sqrt_coeff)}*exp({exp1_coeff_str}*z)))")
print("")

# Print the explanation of the variables
print("Where:")
print("  - `rho(z)` is the gas density at the normalized position z.")
print("  - `rho_0` is the ambient (unshocked) density.")
print("  - `z = x / L` is the position `x` normalized by the characteristic length scale `L`.")
print("  - The length scale is `L = kappa / (rho_0 * M * c_0 * C_v)`.")
print("  - `exp()` is the exponential function, and `sqrt()` is the square root function.")
print("")

# Print the breakdown of the numbers in the equation
print("The specific numbers appearing in this final equation are:")
print(f"  - Numerator: {int(num)}")
print(f"  - Constant term in the denominator: {int(den_c1)}")
print(f"  - Coefficient of z in the exponents: {int(exp1_coeff_num)}/{int(exp1_coeff_den)} and {int(exp2_coeff_num)}/{int(exp2_coeff_den)}")
print(f"  - Coefficient of the second exponential term inside the square root: {int(sqrt_coeff)}")
