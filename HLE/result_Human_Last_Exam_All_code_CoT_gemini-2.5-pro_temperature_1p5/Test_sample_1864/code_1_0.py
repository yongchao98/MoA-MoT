def calculate_flow_controls():
    """
    Calculates the number of essential controls for a flow cytometry experiment.

    This function determines the total number of controls needed by summing up:
    1. The unstained control.
    2. Single-stain controls (one for each channel).
    3. Fluorescence Minus One (FMO) controls (one for each channel).
    """
    
    # Number of channels/fluorochromes in the experiment
    num_channels = 5
    
    # 1. Unstained Control
    unstained_controls = 1
    
    # 2. Single-Stain Compensation Controls (one for each channel)
    single_stain_controls = num_channels
    
    # 3. Fluorescence Minus One (FMO) Controls (one for each channel)
    fmo_controls = num_channels
    
    # Calculate the total number of controls
    total_controls = unstained_controls + single_stain_controls + fmo_controls
    
    # Explain the calculation and print the result
    print("To perform a rigorous 5-color flow cytometry experiment, you need the following controls:")
    print(f"- Unstained Controls: {unstained_controls}")
    print(f"- Single-Stain Compensation Controls: {single_stain_controls}")
    print(f"- Fluorescence Minus One (FMO) Controls: {fmo_controls}")
    print("\nThe total number of essential controls is calculated as follows:")
    print(f"Total = (Unstained) + (Single-Stain) + (FMO)")
    print(f"Equation: {unstained_controls} + {single_stain_controls} + {fmo_controls} = {total_controls}")
    print(f"\nTherefore, you should prepare {total_controls} essential controls for your experiment.")

if __name__ == '__main__':
    calculate_flow_controls()