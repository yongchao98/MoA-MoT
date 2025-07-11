def calculate_steric_virial_coefficient():
    """
    Calculates the second osmotic virial coefficient from steric-only behavior
    for a typical monoclonal antibody.
    """

    # The partial specific volume (v_bar) for a typical monoclonal antibody
    # is approximately 0.73 mL/g. This is a standard value used in biophysical calculations.
    v_bar = 0.73  # units: mL/g

    # The steric-only contribution to the second osmotic virial coefficient (B22)
    # is commonly estimated using the hard-sphere approximation formula: B22_steric = 4 * v_bar.
    # This value represents the contribution of the protein's excluded volume.
    factor = 4
    b22_steric = factor * v_bar

    # Print the explanation and the calculation, showing each number in the equation.
    print("The second osmotic virial coefficient from steric-only behavior (B22, steric) is calculated based on the protein's partial specific volume (v_bar).")
    print("The standard approximation for hard-sphere-like behavior is B22, steric = 4 * v_bar.")
    print("\nUsing a typical v_bar for a monoclonal antibody of {} mL/g:".format(v_bar))
    # Using an f-string to format the output equation clearly.
    print(f"B22, steric = {factor} * {v_bar} = {b22_steric:.3f} mL/g")

calculate_steric_virial_coefficient()