def solve_r0f_expression():
    """
    This function constructs and prints the symbolic expression for R0f.
    """
    # Define the variables as strings for symbolic representation
    b = "b"
    pg = "pg"
    c = "c"
    pt = "pt"
    gamma_t = "𝛾t"
    mu_t = "𝜇t"
    mu_g = "𝜇g"

    # Construct the numerator and denominator of the expression
    numerator = f"{b} * {pg} * {c} * {pt}"
    denominator = f"(({gamma_t} + {mu_t}) * {mu_g})"

    # Format the final expression string
    r0f_expression = f"R0f = {numerator} / {denominator}"

    # Print the final expression
    print(r0f_expression)

# Execute the function to display the answer
solve_r0f_expression()