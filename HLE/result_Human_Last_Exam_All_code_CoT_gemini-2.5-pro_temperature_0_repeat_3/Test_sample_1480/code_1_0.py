def calculate_rule_activation():
    """
    Calculates the rule activation level for a fuzzy set using a t-norm.
    """
    # Given values
    # Phi_k(x') is the firing level of the antecedent part of the rule.
    phi_k_x_prime = 0.7
    # μG_k_j(y_j) is the membership degree of the consequent part of the rule.
    mu_G_k_j_y_j = 0.9

    # A t-norm is used to combine the antecedent and consequent.
    # The most common t-norm is the minimum operator.
    activation_level = min(phi_k_x_prime, mu_G_k_j_y_j)

    # Print the explanation and the final equation with the result
    print("The rule activation level is calculated by applying a t-norm (minimum operator) to the antecedent firing level and the consequent membership degree.")
    print(f"Activation Level = min(Phi_k(x'), μG_k_j(y_j))")
    print(f"Activation Level = min({phi_k_x_prime}, {mu_G_k_j_y_j})")
    print(f"The calculated rule activation level is: {activation_level}")

if __name__ == "__main__":
    calculate_rule_activation()