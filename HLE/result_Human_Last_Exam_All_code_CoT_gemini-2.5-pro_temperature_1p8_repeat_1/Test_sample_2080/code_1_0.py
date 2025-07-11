import sys

def calculate_virial_coefficients():
    """
    This script calculates the steric-only second osmotic virial coefficient for a monoclonal antibody (mAb)
    and breaks down the contributions to the total experimental value.
    """

    # --- Parameters ---
    # The experimental total second osmotic virial coefficient (B22_total) provided in the problem.
    # Unit: mL/g
    b22_total = -7.585

    # The partial specific volume (v_bar) for a typical monoclonal antibody (mAb).
    # This is a standard literature value, assumed here as it was not provided.
    # Unit: mL/g
    v_bar = 0.73

    # --- Calculation ---
    # The steric-only contribution (B22_steric) is calculated from the hard-sphere model.
    # It represents the repulsive force from the molecule's own volume.
    b22_steric = 4 * v_bar

    # For context, we can now determine the electrostatic contribution.
    # It is the difference between the total and the steric components.
    b22_electrostatic = b22_total - b22_steric

    # --- Output ---
    print("--- Analysis of the Second Osmotic Virial Coefficient (B22) ---")
    print(f"The second osmotic virial coefficient from steric-only behavior (B22,steric) is calculated using the hard-sphere approximation.\n")
    print("The formula is:")
    print("  B22,steric = 4 * v_bar (partial specific volume)\n")
    print(f"Using the standard v_bar for a mAb of {v_bar} mL/g, the calculation is:")

    # Output the final equation with each number, as requested.
    print("  B22,steric = 4 * " + f"{v_bar:.3f}" + " mL/g" + f" = {b22_steric:.3f}" + " mL/g")

    print("\nThis positive value represents a repulsive interaction due to molecular volume.")

    # Writing the final answer to a variable for extraction.
    # This part will not be printed to the user console.
    final_answer = f"{b22_steric:.3f}"

    # Adding context using the provided experimental data
    print("\n--- Context: Breakdown of Total B22 ---")
    print(f"  Total B22 (experimental)  : {b22_total:8.3f} mL/g")
    print(f"  Steric B22 (calculated)   : {b22_steric:8.3f} mL/g (Repulsive)")
    print(f"  Electrostatic B22 (derived) : {b22_electrostatic:8.3f} mL/g (Attractive)")
    print("\nThe large negative electrostatic component confirms the strong attractive interactions mentioned in the problem description.")
    
    return final_answer


if __name__ == '__main__':
    calculate_virial_coefficients()
