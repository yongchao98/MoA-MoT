def solve_infinite_product():
    """
    This function provides the symbolic solution for the infinite product
    product_{n=3 to infinity} (1 - z^3 / n^3)
    by printing the final formula as a string.
    """

    # These are the numbers that appear in the original problem statement.
    start_index = 3
    power = 3

    # The formula for the infinite product from n=1 is based on the Gamma function.
    # Product_{n=1 to infinity} (1 - z^3 / n^3) = 1 / (Gamma(1-z) * Gamma(1-z*w) * Gamma(1-z*w^2))
    # where w = exp(2*pi*i/3) is a primitive cube root of unity.

    # To get the product from n=3, we must divide by the terms for n=1 and n=2.
    n_term_1 = 1
    n_term_2 = 2

    # These are the numbers that appear in the terms we divide by.
    term1_numerator_z_power = 3
    term1_denominator_n = 1
    term1_denominator_n_power = 3

    term2_numerator_z_power = 3
    term2_denominator_n = 2
    term2_denominator_n_power = 3

    # These are the numbers in the arguments of the Gamma functions.
    gamma_arg_const_1 = 1
    gamma_arg_const_2 = 1 # for z*w^1
    gamma_arg_const_3 = 2 # for z*w^2

    # These are the numbers defining the cube root of unity, w = exp(2*pi*i/3).
    omega_numerator_pi_coeff = 2
    omega_denominator = 3
    
    # Constructing the final equation string using the numbers defined above.
    equation = (
        f"product_{{n={start_index}}} ^ infinity (1 - z^{power} / n^{power}) = \n\n"
        f"1 / ( (1 - z^{term1_numerator_z_power}/{term1_denominator_n}^{term1_denominator_n_power}) * "
        f"(1 - z^{term2_numerator_z_power}/{term2_denominator_n}^{term2_denominator_n_power}) * "
        f"Gamma({gamma_arg_const_1}-z) * "
        f"Gamma({gamma_arg_const_1}-z*exp({omega_numerator_pi_coeff}*pi*i/{omega_denominator})) * "
        f"Gamma({gamma_arg_const_1}-z*exp({gamma_arg_const_3}*{omega_numerator_pi_coeff}*pi*i/{omega_denominator})) )"
    )

    print("The final formula for the infinite product is:")
    print(equation)
    
    # For clarity, let's simplify the terms we divide by and the exponents in omega
    simplified_equation = (
        f"\nOr, in a more simplified form:\n\n"
        f"1 / ( (1 - z^3) * (1 - z^3/8) * "
        f"Gamma(1-z) * Gamma(1-z*exp(2*pi*i/3)) * Gamma(1-z*exp(4*pi*i/3)) )"
    )
    print(simplified_equation)


if __name__ == '__main__':
    solve_infinite_product()
