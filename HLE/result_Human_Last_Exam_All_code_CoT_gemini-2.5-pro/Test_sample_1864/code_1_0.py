import sys

def calculate_flow_controls():
    """
    Calculates and explains the number of essential controls for a
    five-color flow cytometry sorting experiment using streptavidin beads.
    """
    # Define the number of fluorochromes/channels from the user's query
    num_fluorochromes = 5
    
    # 1. Unstained Control: For setting baseline fluorescence.
    unstained_control = 1
    
    # 2. Single-Stain Compensation Controls: One for each fluorochrome.
    compensation_controls = num_fluorochromes
    
    # 3. Bead-Specific Control: For assessing non-specific binding to streptavidin beads.
    bead_control = 1
    
    # 4. Fluorescence Minus One (FMO) Controls: Crucial for accurate gating in sorting.
    fmo_controls = num_fluorochromes
    
    # Calculate the total number of essential controls
    total_controls = unstained_control + compensation_controls + bead_control + fmo_controls
    
    # --- Output Section ---
    print("To ensure high-quality data from your five-color sorting experiment, you should prepare the following essential controls:")
    print(f"\n1. Unstained Control: {unstained_control}")
    print("   - Purpose: To measure autofluorescence and set baseline voltages.")
    
    print(f"\n2. Single-Stain (Compensation) Controls: {compensation_controls}")
    print("   - Purpose: To calculate the compensation matrix and correct for spectral spillover from each of the five fluorochromes.")

    print(f"\n3. Bead Non-Specific Binding Control: {bead_control}")
    print("   - Purpose: To measure background signal from fluorescent reagents binding directly to the streptavidin beads.")

    print(f"\n4. Fluorescence Minus One (FMO) Controls: {fmo_controls}")
    print("   - Purpose: To accurately set gates for identifying positive populations, which is critical for sort purity.")
    
    print("\n--------------------------------------------------")
    print("Total number of essential controls is the sum of these individual controls:")
    print(f"{unstained_control} (Unstained) + {compensation_controls} (Compensation) + {bead_control} (Bead) + {fmo_controls} (FMO) = {total_controls}")
    
    # Final answer in the specified format
    # The 'file=sys.stdout' ensures it prints to the standard output, which is the default.
    print(f"\n<<<{total_controls}>>>", file=sys.stdout)

if __name__ == '__main__':
    calculate_flow_controls()