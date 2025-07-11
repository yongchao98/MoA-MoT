def print_inner_product_expression():
    """
    This function prints the mathematical expression for the inner product (phi, D_phi)
    in the functional integral representation of the partition function for a neutral
    scalar field at finite temperature.

    Symbols:
    - phi: The neutral scalar field phi(tau, x).
    - D: The Euclidean Klein-Gordon operator.
    - (f, g): The inner product, defined as an integral over Euclidean spacetime.
    - tau: Imaginary time, integrated from 0 to beta (inverse temperature).
    - d^d x: Integral over d spatial dimensions.
    - partial_tau: The partial derivative with respect to tau.
    - nabla: The spatial gradient operator.
    - m: The mass of the scalar field.
    """

    # Define the components of the equation as strings
    inner_product_lhs = "(phi, D_phi)"
    integral_part = "∫_0^β dτ ∫ d^d x"
    
    # The integrand contains the kinetic and mass terms
    kinetic_term_time = "(∂_τ φ)²"
    kinetic_term_space = "(∇φ)²"
    mass_term = "m²φ²"
    
    integrand = f"[{kinetic_term_time} + {kinetic_term_space} + {mass_term}]"

    # Construct the final equation string
    final_equation = f"{inner_product_lhs} = {integral_part} {integrand}"

    # Print the equation
    print("In finite-temperature field theory for a neutral scalar field, the inner product (φ, D_φ) is given by:")
    print(final_equation)
    
    # Print the definition of the operator D
    print("\nWhere the operator D is the Euclidean Klein-Gordon operator:")
    print("D = -∂_τ² - ∇² + m²")


if __name__ == "__main__":
    print_inner_product_expression()