import math

def generate_vertical_slice_formulation():
    """
    Calculates and displays the mathematical formulation for a vertical
    cross-section of an Interval Type-3 Membership Function (IT3 MF).

    This formulation uses a Gaussian function to represent the upper bound
    of the membership grade for the secondary variable 'u', given a fixed
    primary input 'x_prime'.
    """

    # 1. Define a fixed value for the primary input variable 'x'.
    x_prime = 0.5

    # 2. Define the parameters of the Gaussian function (height, center, sigma)
    #    as functions of the primary input 'x'. These are examples to illustrate
    #    how the slice shape changes with 'x'.
    # For the center of the Gaussian: c(x) = 0.4*x + 0.3
    # For the standard deviation (width) of the Gaussian: sigma(x) = -0.2*x + 0.4
    # For the height of the Gaussian: h(x) = 0.95 (constant in this example)

    center_val = 0.4 * x_prime + 0.3
    sigma_val = -0.2 * x_prime + 0.4
    height_val = 0.95

    # 3. Explain the context and print the calculated parameter values.
    print(f"For a fixed primary input x = {x_prime}, the vertical cross-section is a Type-2 Fuzzy Set.")
    print("Its upper bound can be modeled by a Gaussian function μ(u) = h * exp(-((u - c) / σ)^2).")
    print("\nCalculated parameters for this cross-section:")
    
    # Print each number that will be in the final equation
    print(f"Height (h): {height_val}")
    print(f"Center (c): {center_val}")
    print(f"Standard Deviation (σ): {sigma_val}")
    
    # 4. Construct and print the final mathematical formulation.
    # The 'u' remains a variable in the final formula.
    # Note: The '2' in the exponent is also a number in the equation.
    equation_exponent = 2
    
    equation_str = (
        f"μ(u) = {height_val:.2f} * exp(-((u - {center_val:.2f}) / {sigma_val:.2f})^ {equation_exponent})"
    )

    print("\nThe final mathematical formulation for the upper bound is:")
    print(equation_str)

    # Return the final answer in the required format.
    # The final answer is the mathematical formulation string.
    final_answer = f"<<<{equation_str}>>>"
    print(final_answer)

generate_vertical_slice_formulation()
