def calculate_activation_level():
    """
    Calculates the rule activation level for a fuzzy set using a t-norm.
    """
    # Given membership grades from the problem description
    phi_k_x_prime = 0.7
    mu_G_k_j_y_j = 0.9

    # A common t-norm is the minimum operator. The rule activation level is the
    # t-norm of the antecedent membership grades.
    activation_level = min(phi_k_x_prime, mu_G_k_j_y_j)

    # Output the explanation and the equation as requested
    print("The rule activation level is calculated by applying a t-norm (e.g., minimum) to the antecedent membership grades.")
    print(f"Activation Level = min({phi_k_x_prime}, {mu_G_k_j_y_j})")
    print(f"Result: {activation_level}")

calculate_activation_level()
<<<0.7>>>