def solve_rule_activation():
    """
    Calculates the rule activation level for a fuzzy set using a t-norm.

    The problem provides the firing strength from the antecedent (Phi_k) and
    the membership grade from the consequent (mu_G). The rule activation level
    is found by applying a t-norm to these two values. We will use the
    most common t-norm, the minimum operator.
    """
    # Given values
    phi_k = 0.7
    mu_G = 0.9

    # The t-norm operation (minimum)
    # Activation Level = min(Phi_k(x'), μG_k_j(y_j))
    activation_level = min(phi_k, mu_G)

    # Print the final equation with the numbers and the result
    print(f"min({phi_k}, {mu_G}) = {activation_level}")

solve_rule_activation()
<<<0.7>>>