import sympy as sp

def evaluate_1s_integral():
    """
    Evaluates the integral <phi_1s| 1/r |phi_1s> for a Slater-Type Orbital
    using symbolic mathematics with sympy.
    """
    # Step 1: Define the symbolic variables
    # r: radial distance
    # theta: polar angle
    # phi: azimuthal angle
    # zeta: orbital exponent (a positive real number)
    r, theta, phi = sp.symbols('r theta phi', real=True)
    zeta = sp.symbols('zeta', positive=True)

    print("Evaluating the integral <phi_1s| 1/r |phi_1s> for a 1s Slater orbital.")
    print("The orbital is phi_1s(r) = N * exp(-zeta*r). The operator is 1/r.")
    print("-" * 60)

    # Step 2: Calculate the components of the integral expression.
    # The full integral is: Integral( (N * exp(-zeta*r))^2 * (1/r) * r^2 * sin(theta) dr dtheta dphi )
    # This can be broken down into three parts:
    # Part 1: The squared normalization constant, N^2
    # Part 2: The radial integral, Integral( r * exp(-2*zeta*r) dr )
    # Part 3: The angular integral, Integral( sin(theta) dtheta dphi )

    # Part 1: Calculate N^2. N is found from the normalization condition:
    # N^2 * Integral( (exp(-zeta*r))^2 * r^2 * sin(theta) dr dtheta dphi ) = 1
    unnormalized_psi = sp.exp(-zeta * r)
    norm_integral = sp.integrate(
        unnormalized_psi**2 * r**2 * sp.sin(theta),
        (r, 0, sp.oo),
        (theta, 0, sp.pi),
        (phi, 0, 2 * sp.pi)
    )
    N_squared = 1 / norm_integral

    # Part 2: The radial part of the expectation value integral.
    # After combining terms, the radial part is Integral( r * exp(-2*zeta*r) dr, {r, 0, oo} )
    radial_integral = sp.integrate(r * sp.exp(-2 * zeta * r), (r, 0, sp.oo))

    # Part 3: The angular part of the expectation value integral.
    # This is Integral( sin(theta) dtheta dphi ) over the sphere.
    angular_integral = sp.integrate(
        sp.sin(theta),
        (theta, 0, sp.pi),
        (phi, 0, 2 * sp.pi)
    )

    # Step 3: Combine the parts to get the final equation.
    # Result = N^2 * radial_integral * angular_integral
    final_result = N_squared * radial_integral * angular_integral
    
    # Step 4: Print the final equation with each component shown.
    print("The expectation value is the product of three terms:")
    print(f"1. Squared Normalization Constant (N^2): {N_squared}")
    print(f"2. Radial Integral: {radial_integral}")
    print(f"3. Angular Integral: {angular_integral}")
    print("-" * 60)
    print("The final equation is the product of these terms:")
    # The instruction is to output each number in the final equation.
    # Here, we show the symbolic expression for each term.
    print(f"Result = ({N_squared}) * ({radial_integral}) * ({angular_integral})")

    # Print the simplified final result
    print("\nAfter simplification, the final result is:")
    print(f"Result = {final_result}")

if __name__ == '__main__':
    evaluate_1s_integral()
