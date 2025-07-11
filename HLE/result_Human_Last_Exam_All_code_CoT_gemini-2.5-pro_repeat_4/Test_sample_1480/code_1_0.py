def calculate_rule_activation():
    """
    Calculates the rule activation level for a fuzzy rule using a t-norm.

    The rule activation level is found by combining the antecedent's firing strength
    with the consequent's membership value. This is typically done using a t-norm
    operation, which models a fuzzy 'AND'. The most common t-norm is the minimum
    operator, which we will use here.
    """
    # Given values
    phi_k = 0.7      # Firing level of the antecedent: Phi_k(x')
    mu_G_k_j = 0.9   # Membership value in the consequent set: μG_k_j(y_j)

    # Using the minimum t-norm for the calculation
    # Activation Level = t-norm(phi_k, mu_G_k_j) = min(phi_k, mu_G_k_j)
    activation_level = min(phi_k, mu_G_k_j)

    # Output the result in a clear format, showing the equation
    print("The rule activation level is calculated using the minimum t-norm.")
    print(f"Activation Level = min({phi_k}, {mu_G_k_j})")
    print(f"Activation Level = {activation_level}")

# Execute the function
calculate_rule_activation()