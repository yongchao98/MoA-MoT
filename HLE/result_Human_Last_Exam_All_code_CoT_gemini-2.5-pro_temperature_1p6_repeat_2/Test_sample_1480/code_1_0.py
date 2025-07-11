def calculate_activation_level():
    """
    Calculates the rule activation level for an interval type-3 fuzzy set
    using a minimum t-norm.
    """
    # Given values
    phi_k = 0.7
    mu_G_k_j = 0.9

    # The rule activation level is found by applying a t-norm to the antecedent
    # firing level and the consequent membership degree.
    # We will use the minimum t-norm: activation_level = min(phi_k, mu_G_k_j)
    activation_level = min(phi_k, mu_G_k_j)

    print("The rule activation level is calculated using a t-norm (e.g., minimum) on the antecedent firing level and the consequent membership degree.")
    print(f"Antecedent Firing Level (Phi_k(x′)): {phi_k}")
    print(f"Consequent Membership Degree (μG_k_j(y_j)): {mu_G_k_j}")
    print("\nCalculating the rule activation level using the minimum t-norm:")
    print(f"Activation Level = min({phi_k}, {mu_G_k_j})")
    print(f"Final Rule Activation Level = {activation_level}")

calculate_activation_level()