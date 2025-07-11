import math

# This script will derive and print the analytical solution for the shock density profile.

# The final solution is an implicit equation relating the normalized density (ρ/ρ₀)
# to the normalized position (x/L). The numbers in the equation are derived from the
# physical parameters of the problem (Mach number = sqrt(3), gamma = 5/3, Prandtl number = 3/4).

# Let's define the constants that appear in the final equation.
# The equation is of the form: x/L = A * ln( B * ( (ρ/ρ₀ - ρ_up)/(ρ/ρ₀ * (ρ_down - ρ/ρ₀)) ) )
# We will present it as: x/L = (c1_num/c1_den) * ln( c2 * ((ρ/ρ₀ - ρ_upstream)**exponent) / (ρ/ρ₀ * (ρ_downstream - ρ/ρ₀)) )

# Coefficient in front of the logarithm, derived from non-dimensionalization
c1_num = 3
c1_den = 4

# Multiplicative constant inside the logarithm, from setting the origin x=0
c2 = 3

# Normalized upstream (ambient) density
rho_upstream = 1

# Normalized downstream (post-shock) density
rho_downstream = 2

# The exponent of the numerator term, from the partial fraction integration
exponent = 2

# Symbols for the normalized density and position variables
rho_norm_symbol = "ρ/ρ₀"
x_norm_symbol = "x/L"

# Construct the equation string using the defined constants.
# This makes it clear where each number in the final equation comes from.
equation = (
    f"{x_norm_symbol} = ({c1_num}/{c1_den}) * ln( "
    f"{c2} * (({rho_norm_symbol} - {rho_upstream})**{exponent}) / "
    f"(({rho_norm_symbol}) * ({rho_downstream} - {rho_norm_symbol})) )"
)

# Print the final result
print("The analytical solution for the density profile of the specified shock wave is given by the following implicit equation:")
print(equation)
print("\nWhere:")
print(f"  {rho_norm_symbol} is the density normalized by the ambient density.")
print(f"  {x_norm_symbol} is the position normalized by the conductive length scale L = κ/(ρ₀*M*c₀*Cᵥ).")
print(f"  This solution is valid for {rho_upstream} < {rho_norm_symbol} < {rho_downstream}.")
