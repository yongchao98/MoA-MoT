def calculate_flow_controls():
    """
    Calculates the number of essential controls for a multi-color flow cytometry experiment.
    """
    # Number of fluorescent channels/colors in the experiment
    num_colors = 5

    # 1. Unstained Control: One sample with just the beads.
    unstained_controls = 1
    print(f"Unstained Controls: {unstained_controls}")

    # 2. Single-Stain Controls: One for each color to set compensation.
    single_stain_controls = num_colors
    print(f"Single-Stain Compensation Controls: {single_stain_controls}")

    # 3. Fluorescence Minus One (FMO) Controls: One for each color to set accurate gates.
    fmo_controls = num_colors
    print(f"Fluorescence Minus One (FMO) Controls: {fmo_controls}")

    # Calculate the total number of controls
    total_controls = unstained_controls + single_stain_controls + fmo_controls

    # Print the final calculation and result
    print("\n--- Total Calculation ---")
    print(f"Total Essential Controls = {unstained_controls} (Unstained) + {single_stain_controls} (Single-Stains) + {fmo_controls} (FMOs)")
    print(f"Total Essential Controls = {total_controls}")

if __name__ == "__main__":
    calculate_flow_controls()