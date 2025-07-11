import math

# The problem asks for the analytical solution for the density profile.
# This script will print the final equation in a readable format.
# All numbers in the equation are derived from the problem's specific
# parameters (Mach number = sqrt(3), gamma = 5/3, Prandtl number = 3/4).

# Constants in the final equation
val_1_A = 1
val_1_B = 1
val_1_C = 1
val_3 = 3
val_2 = 2
val_3_exp = 3
val_0_5 = 0.5

print("The analytical solution for the density profile is given by the equation:")
print("")
# Using an f-string to embed the numbers directly into the equation string.
print(f"rho/rho_0 = {val_1_A} + {val_1_B} / ( {val_1_C} + {val_3} * exp(-({val_2}/{val_3_exp}) * x/L) )**{val_0_5}")
print("")
print("where:")
print("  rho       is the local gas density.")
print("  rho_0     is the ambient density upstream of the shock.")
print("  x         is the position along the direction of propagation.")
print("  L         is the ambient conductive length scale, defined as:")
print("            L = kappa / (rho_0 * M * c_0 * C_v)")
print("  kappa     is the thermal conductivity.")
print(f"  M         is the Mach number, M = sqrt(3) ~ {math.sqrt(3):.3f}")
print("  c_0       is the ambient speed of sound.")
print("  C_v       is the specific heat at constant volume.")
print("  exp()     is the exponential function.")
