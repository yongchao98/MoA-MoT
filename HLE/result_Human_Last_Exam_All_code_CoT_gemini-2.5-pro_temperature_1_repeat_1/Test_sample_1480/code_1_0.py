def calculate_rule_activation_level():
    """
    Calculates the rule activation level for an interval type-3 fuzzy set
    using a minimum t-norm operation.
    """
    # Given values for the firing strengths
    phi_k_x_prime = 0.7
    mu_G_k_j_y_j = 0.9

    # The rule activation level is calculated by applying a t-norm to the
    # antecedent and consequent firing strengths. We will use the minimum t-norm,
    # which is the most common choice.
    # Activation Level = min(Phi_k(x'), μG_k_j(y_j))
    activation_level = min(phi_k_x_prime, mu_G_k_j_y_j)

    # Print the details of the calculation
    print("Calculating the rule activation level using the minimum t-norm:")
    print(f"Rule Activation Level = min(Phi_k(x'), μG_k_j(y_j))")
    print(f"Rule Activation Level = min({phi_k_x_prime}, {mu_G_k_j_y_j})")
    print(f"The final rule activation level is: {activation_level}")

# Execute the function
calculate_rule_activation_level()