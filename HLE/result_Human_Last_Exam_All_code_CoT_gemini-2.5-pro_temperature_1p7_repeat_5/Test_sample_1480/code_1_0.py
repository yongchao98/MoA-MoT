def calculate_activation_level():
    """
    Calculates the rule activation level for a fuzzy set using a minimum t-norm.
    """
    # Given values for the problem
    phi_k = 0.7
    mu_G_k_j = 0.9

    # The rule activation level is determined by applying a t-norm operation
    # to the firing strengths. We will use the most common t-norm, the minimum operator.
    activation_level = min(phi_k, mu_G_k_j)

    # Print the process and the result
    print("Calculating the rule activation level using the minimum t-norm:")
    print(f"Activation Level = min(Phi_k(x′), μG_k_j(y_j))")
    print(f"Activation Level = min({phi_k}, {mu_G_k_j})")
    print(f"Result: {activation_level}")

calculate_activation_level()