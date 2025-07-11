def print_r0f_expression():
    """
    This function constructs and prints the symbolic expression for R0f.
    It follows the user's request to output each variable in the final equation.
    """
    # Define the variables as strings for symbolic representation
    b = "b"
    c = "c"
    pg = "pg"
    pt = "pt"
    gamma_t = "𝛾t"
    mu_t = "𝜇t"
    mu_g = "𝜇g"

    # Print the final equation piece by piece to show all components
    print("The final expression for R0f is:")
    print("R0f", "=", end=" ")
    print("(", b, "*", c, "*", pg, "*", pt, ")", end=" ")
    print("/", end=" ")
    print("(", "(", gamma_t, "+", mu_t, ")", "*", mu_g, ")")

# Execute the function to print the expression
print_r0f_expression()