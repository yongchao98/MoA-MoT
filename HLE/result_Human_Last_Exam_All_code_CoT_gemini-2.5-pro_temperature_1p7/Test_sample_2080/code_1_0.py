import math

def calculate_steric_virial_coefficient():
    """
    Calculates the second osmotic virial coefficient from steric-only behavior
    for a monoclonal antibody.
    """
    # The second osmotic virial coefficient from steric-only behavior (B_2,steric)
    # is calculated using the hard-sphere model.
    # The formula is B_2,steric = 4 * v_sp, where v_sp is the partial specific volume.

    # A typical partial specific volume (v_sp) for a monoclonal antibody is 0.73 mL/g.
    partial_specific_volume = 0.73  # in mL/g

    # The hard-sphere model factor is 4.
    hard_sphere_factor = 4

    # The total measured B22 value is given as -7.585 mL/g.
    # This value represents the sum of steric and all other interactive forces
    # and is provided for context, not for the calculation of the steric-only component.
    b22_total = -7.585 # in mL/g

    # Calculate the steric-only component
    b2_steric = hard_sphere_factor * partial_specific_volume

    print("The second osmotic virial coefficient from Steric-only behavior (B_2,steric) is calculated from the hard-sphere model.")
    print(f"A typical partial specific volume (v_sp) for a monoclonal antibody is {partial_specific_volume} mL/g.")
    print("\nThe calculation is based on the formula: B_2,steric = factor * v_sp")
    print(f"\nFinal Equation: {hard_sphere_factor} * {partial_specific_volume} = {b2_steric:.3f}")
    print(f"The resulting steric-only virial coefficient is {b2_steric:.3f} mL/g.")

calculate_steric_virial_coefficient()