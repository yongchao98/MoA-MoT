def calculate_activation_level():
    """
    Calculates the rule activation level for an interval type-3 fuzzy set
    using a minimum t-norm operation.
    """
    # Given values
    phi_k = 0.7
    mu_G_k_j = 0.9

    # The rule activation level (w) is calculated by applying a t-norm
    # to the antecedent's firing level and the consequent's membership grade.
    # We will use the minimum t-norm, which is the most common choice.
    # Formula: w = min(Phi_k(x'), μG_k_j(y_j))
    activation_level = min(phi_k, mu_G_k_j)

    # Print the explanation and the final equation with the result
    print("The rule activation level (w) is calculated using the minimum t-norm:")
    print(f"w = min(Phi_k(x'), μG_k_j(y_j))")
    print(f"w = min({phi_k}, {mu_G_k_j})")
    print(f"The final rule activation level is: {activation_level}")

calculate_activation_level()
<<<0.7>>>