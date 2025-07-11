def calculate_activation_level():
    """
    Calculates the rule activation level for an interval type-3 fuzzy set
    using a minimum t-norm operation.
    """
    # Given values
    phi_k = 0.7  # Antecedent firing level
    mu_g = 0.9   # Consequent membership degree

    # The rule activation level is calculated using a t-norm.
    # The minimum operator is a common t-norm.
    # Activation Level = t-norm(phi_k, mu_g)
    activation_level = min(phi_k, mu_g)

    # Print the equation and the result
    print("The rule activation level is determined by applying a t-norm (minimum) to the antecedent firing level and the consequent membership degree.")
    print(f"Rule Activation Level = min(Phi_k(x'), μG_k_j(y_j))")
    print(f"Rule Activation Level = min({phi_k}, {mu_g}) = {activation_level}")

# Execute the function
calculate_activation_level()