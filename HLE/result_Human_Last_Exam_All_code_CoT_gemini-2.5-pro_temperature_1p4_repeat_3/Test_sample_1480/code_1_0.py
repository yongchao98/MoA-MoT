import numpy as np

def calculate_activation_level():
    """
    Calculates the rule activation level using a t-norm operation.

    The problem provides:
    - Antecedent firing level, Phi_k(x'): 0.7
    - Consequent membership grade, muG_k_j(y_j): 0.9

    The t-norm operation combines these values. We will use the most
    common t-norm, which is the minimum function.
    """
    
    # Step 1: Define the given values
    antecedent_firing_level = 0.7
    consequent_membership_grade = 0.9

    # Step 2: Apply the minimum t-norm to calculate the activation level
    activation_level = min(antecedent_firing_level, consequent_membership_grade)

    # Step 3: Print the final equation and the result
    print("The rule activation level is calculated using the minimum t-norm as follows:")
    print(f"Activation Level = min({antecedent_firing_level}, {consequent_membership_grade}) = {activation_level}")

calculate_activation_level()