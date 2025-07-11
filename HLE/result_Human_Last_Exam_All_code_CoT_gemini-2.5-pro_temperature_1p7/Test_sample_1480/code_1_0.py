import sys

def calculate_rule_activation():
    """
    Calculates the rule activation level for an interval type-3 fuzzy set
    using a minimum t-norm operation.
    """
    # Given membership values
    phi_k = 0.7
    mu_G = 0.9

    # The rule activation level is determined by applying a t-norm to the
    # membership values in the rule's antecedent. The most common t-norm
    # is the minimum (min) function.
    activation_level = min(phi_k, mu_G)

    # Print the equation and the final result.
    # We use sys.stdout.write to ensure there's no extra formatting.
    sys.stdout.write("The rule activation level is calculated using the minimum t-norm:\n")
    sys.stdout.write(f"Activation Level = min({phi_k}, {mu_G})\n")
    sys.stdout.write(f"Activation Level = {activation_level}\n")

calculate_rule_activation()