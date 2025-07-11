import math

# This script prints the analytical solution for the specified shock wave density profile.

# The derived relationship involves the normalized density `r` and the normalized position `x_tilde`.
# r = rho / rho_0, where rho is the local density and rho_0 is the ambient density.
# x_tilde = x / L, where L is the characteristic length scale.

# The problem is for a Mach sqrt(3) shock in a monatomic ideal gas with Pr = 3/4.
# This gives upstream and downstream density ratios of r_0 = 1 and r_1 = 2.

# The analytical solution for the profile is given by the following equation,
# centered at r = 1.5.

# Constants in the final equation
c1_num = 3
c1_den = 4
c2 = 3
c3 = 1
c4 = 2

print("The analytical solution relating the normalized density r = rho/rho_0")
print("and the normalized position x_tilde = x/L is:")
print("")
print(f"x_tilde = ({c1_num}/{c1_den}) * ln[ {c2} * (r - {c3})^2 / (r * ({c4} - r)) ]")
print("")
print("Where:")
print("  r         is the density in units of the ambient density (rho_0).")
print("            It varies from 1 (far upstream) to 2 (far downstream).")
print("  x_tilde   is the position in units of the ambient conductive length scale,")
print("            L = kappa/(rho_0 * M * c_0 * C_v).")
print("  ln        is the natural logarithm.")
print("  kappa     is the thermal conductivity.")
print("  rho_0     is the ambient density.")
print("  M         is the shock Mach number (sqrt(3)).")
print("  c_0       is the ambient sound speed.")
print("  C_v       is the specific heat at constant volume.")
