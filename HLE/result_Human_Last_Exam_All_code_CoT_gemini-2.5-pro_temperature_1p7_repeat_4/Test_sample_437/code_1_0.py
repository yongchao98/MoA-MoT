import math

def evaluate_electron_nucleus_attraction(zeta):
    """
    Evaluates the electron-nucleus attraction integral <phi_1s| 1/r |phi_1s>
    for a 1s Slater-type orbital (STO) with a given exponent zeta.

    The final result of this integral is simply zeta. This function
    demonstrates the calculation by evaluating each part of the integral's
    analytical solution: (zeta^3/pi) * (4*pi) * (1/(4*zeta^2)).
    """
    if zeta <= 0:
        print("Zeta must be a positive number.")
        return

    print(f"Evaluating the integral for zeta = {zeta:.4f}")

    # Part 1: The square of the normalization constant from phi_1s^2
    norm_squared = zeta**3 / math.pi

    # Part 2: The result from integrating over the angular coordinates (theta, phi)
    angular_integral = 4 * math.pi

    # Part 3: The result from integrating over the radial coordinate (r)
    radial_integral = 1 / (2 * zeta)**2

    # The final result is the product of these three parts
    final_value = norm_squared * angular_integral * radial_integral

    # As requested, printing each number in the final equation
    print("\nThe integral is calculated as the product of three terms:")
    print("Term 1 (from normalization): zeta^3 / pi")
    print("Term 2 (from angular integration): 4 * pi")
    print("Term 3 (from radial integration): 1 / (2 * zeta)^2")
    
    print("\nThe final equation with numerical values is:")
    print(f"({norm_squared:.4f}) * ({angular_integral:.4f}) * ({radial_integral:.4f}) = {final_value:.4f}")

    print("\nVerifying the simplified formula: Result = zeta")
    print(f"The calculated value is indeed {final_value:.4f}, which matches the expected value of zeta = {zeta:.4f}.")


# --- Main execution ---
# Set the value for the orbital exponent zeta.
# For a hydrogen atom, zeta = 1.0.
# For a helium atom, a common approximation is zeta = 1.6875.
# Let's use zeta for hydrogen as an example.
zeta_value = 1.0
evaluate_electron_nucleus_attraction(zeta_value)
