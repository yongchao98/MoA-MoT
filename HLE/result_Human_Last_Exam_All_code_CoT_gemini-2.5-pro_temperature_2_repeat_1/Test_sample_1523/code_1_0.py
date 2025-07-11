def print_inner_product_expression():
    """
    This function constructs and prints the mathematical expression for the
    inner product (ϕ, D_ϕ) for a neutral scalar field at finite temperature.
    It explicitly details the numerical coefficients as requested.
    """

    # Define the various mathematical components as strings
    inner_product_symbol = "(ϕ, D_ϕ)"
    integral = "∫_0^β dτ ∫ d³x"
    field = "ϕ(τ, x)"

    # Define the terms of the operator D = -∂_τ² - ∇² + m²
    op_term1 = "∂_τ²"
    op_term2 = "∇²"
    op_term3 = "m²"

    # Define the numerical coefficients for each term in the operator
    coeff1 = -1
    coeff2 = -1
    coeff3 = 1

    # Construct the full expression string
    final_equation = (f"{inner_product_symbol} = {integral} {field} * "
                      f"[ ({coeff1})*{op_term1} + ({coeff2})*{op_term2} + ({coeff3})*{op_term3} ] "
                      f"* {field}")

    # Print the final result
    print("In the functional integral representation for a neutral scalar field, the inner product (ϕ, D_ϕ) that appears in the quadratic part of the action is:")
    print(final_equation)
    print("\nThis expression shows the integral of the field ϕ multiplied by the result of the kinetic operator D acting on ϕ.")

    # Explicitly print each number in the equation as requested
    print("\nThe numerical coefficients in the operator D are:")
    print(f"1. The coefficient of the second time-derivative term ({op_term1}) is: {coeff1}")
    print(f"2. The coefficient of the Laplacian operator ({op_term2}) is: {coeff2}")
    print(f"3. The coefficient of the mass-squared term ({op_term3}) is: {coeff3}")


# Execute the function to display the result
print_inner_product_expression()