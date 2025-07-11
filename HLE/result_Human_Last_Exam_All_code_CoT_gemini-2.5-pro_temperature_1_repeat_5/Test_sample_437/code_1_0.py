import math

def evaluate_integral():
    """
    This function evaluates the integral <phi_i| 1/r |phi_j> for 1s Slater-type
    orbitals using the derived analytical formula and prints the calculation.
    """
    # We will evaluate the integral for a general case with two different
    # orbital exponents, zeta_i and zeta_j.
    # As an example, let's choose values that could be used in a variational
    # calculation for a two-electron system.
    zeta_i = 1.0
    zeta_j = 1.2

    # The analytical formula for the integral is:
    # I = 4 * (zeta_i * zeta_j)^(3/2) / (zeta_i + zeta_j)^2
    
    # Calculate the result using the formula
    result = 4 * (zeta_i * zeta_j)**1.5 / (zeta_i + zeta_j)**2

    # Print the final equation with the numbers substituted in, as requested.
    print(f"For 1s Slater orbitals with exponents zeta_i = {zeta_i} and zeta_j = {zeta_j}, the integral <phi_i| 1/r |phi_j> is evaluated as:")
    print(f"4 * ({zeta_i} * {zeta_j})**1.5 / ({zeta_i} + {zeta_j})**2 = {result}")

    # For verification, we also demonstrate the important special case where 
    # zeta_i = zeta_j = zeta. For zeta = 1.0 (e.g., a hydrogen atom), 
    # the integral should simply be equal to zeta.
    zeta = 1.0
    special_result = 4 * (zeta * zeta)**1.5 / (zeta + zeta)**2
    print("\nFor the special case where zeta_i = zeta_j = 1.0, the integral becomes:")
    print(f"4 * ({zeta} * {zeta})**1.5 / ({zeta} + {zeta})**2 = {special_result}")

# Execute the function
evaluate_integral()