def calculate_rule_activation():
    """
    Calculates the rule activation level using a minimum t-norm.
    """
    # Given values
    phi_k = 0.7
    mu_G = 0.9

    # The t-norm operation is used to combine the antecedent firing strength
    # with the consequent membership function. We will use the most common
    # t-norm, which is the minimum operator.
    activation_level = min(phi_k, mu_G)

    # Print the equation and the result
    print(f"The rule activation level is calculated using the minimum t-norm:")
    print(f"min({phi_k}, {mu_G}) = {activation_level}")

calculate_rule_activation()