def solve_and_print_product_formula():
    """
    This function prints the final symbolic formula for the infinite product.
    The formula is derived from the Weierstrass product for the Gamma function.
    The output is formatted to clearly show all the numbers involved in the equation.
    """

    # The simplified formula for the product is:
    # 8 / ( (1 - z^3) * (8 - z^3) * Gamma(1-z) * Gamma(1 - z*exp(2*pi*I/3)) * Gamma(1 - z*exp(4*pi*I/3)) )
    
    # Define the numbers present in the formula to satisfy the prompt's requirement.
    numerator = 8
    term1_const = 1
    term1_power = 3
    term2_const = 8
    term2_power = 3
    gamma_arg_const = 1
    omega_factor_num = 2
    omega_factor_den = 3
    omega_sq_factor_num = 4
    omega_sq_factor_den = 3

    # Print the equation in a formatted, readable way.
    # This print statement fulfills the requirement to "output each number in the final equation".
    print(
        f"  {numerator}\n"
        "------------------------------------------------------------------------------------------------------------------\n"
        f"({term1_const} - z^{term1_power})({term2_const} - z^{term2_power}) * Gamma({gamma_arg_const} - z) * Gamma({gamma_arg_const} - z*exp({omega_factor_num}*pi*I/{omega_factor_den})) * Gamma({gamma_arg_const} - z*exp({omega_sq_factor_num}*pi*I/{omega_sq_factor_den}))"
    )

if __name__ == "__main__":
    solve_and_print_product_formula()