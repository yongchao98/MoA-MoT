def print_infinite_product_solution():
    """
    This function prints the closed-form symbolic solution for the given infinite product.
    The solution is expressed in terms of the Gamma function.
    """
    # Define the numbers that appear in the final equation to fulfill the user request.
    start_index = 3
    power = 3
    denominator_n1_cubed = 1**3
    denominator_n2_cubed = 2**3
    gamma_arg_const = 1
    exp_arg1_num = 2
    exp_arg1_den = 3
    exp_arg2_num = 4
    exp_arg2_den = 3

    # Construct the final expression as a string.
    # Note: 'z' is a complex variable, 'i' is the imaginary unit, 'pi' is the constant pi,
    # 'exp' is the exponential function, and 'Gamma' represents the Gamma function.
    # The denominator terms for n=1 and n=2 are simplified for readability.
    final_expression = (
        f"1 / ( (1 - z**{power}) * (1 - z**{power}/{denominator_n2_cubed}) * "
        f"Gamma({gamma_arg_const} - z) * "
        f"Gamma({gamma_arg_const} - z * exp(i*{exp_arg1_num}*pi/{exp_arg1_den})) * "
        f"Gamma({gamma_arg_const} - z * exp(i*{exp_arg2_num}*pi/{exp_arg2_den})) )"
    )

    print(f"The expression for the infinite product is:")
    print(final_expression)

print_infinite_product_solution()