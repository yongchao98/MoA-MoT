import math

def print_invariant_density():
    """
    This function derives and prints the normalised density of the invariant measure
    for the map T(x) = 1/sqrt(x) mod 1.

    The derivation is based on a change of variable u = sqrt(x), which transforms
    the map into the Gauss map T(u) = 1/u mod 1. The known invariant density
    for the Gauss map is then transformed back to the original variable x.
    """

    # The invariant density for the Gauss map T(u) = 1/u mod 1 is rho_u(u) = 1 / (ln(2)*(1+u)).
    # The transformation of densities from u to x is given by rho_x(x) = rho_u(sqrt(x)) * |du/dx|.
    # Since u = sqrt(x), du/dx = 1 / (2*sqrt(x)).
    # So, rho_x(x) = (1 / (ln(2)*(1+sqrt(x)))) * (1 / (2*sqrt(x))).
    # This simplifies to rho_x(x) = 1 / (2*ln(2)*(sqrt(x) + x)).

    # The final equation for the density rho(x) is:
    # rho(x) = 1 / (C * (x**0.5 + x))
    # where C = 2 * ln(2)

    # Let's define and print the components of the equation.
    numerator = 1
    constant_factor_2 = 2
    constant_factor_ln2 = math.log(2)
    
    print("The normalised density of the invariant measure is rho(x).")
    print("The formula for the density is derived to be:")
    print(f"rho(x) = {numerator} / ({constant_factor_2} * ln(2) * (x^0.5 + x))\n")

    print("As requested, here are the numbers in the final equation:")
    print(f"Numerator: {numerator}")
    
    print("\nDenominator components:")
    print(f"  - Constant factor: {constant_factor_2}")
    print(f"  - Logarithmic constant: ln(2) ≈ {constant_factor_ln2:.15f}")
    print(f"  - The full constant C = 2 * ln(2) ≈ {constant_factor_2 * constant_factor_ln2:.15f}")
    print("  - Variable part: (x^0.5 + x)")

print_invariant_density()