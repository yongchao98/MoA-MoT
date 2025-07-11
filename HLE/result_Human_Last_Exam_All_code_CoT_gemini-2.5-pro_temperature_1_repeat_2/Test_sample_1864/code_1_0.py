def calculate_flow_controls():
    """
    Calculates the number of essential controls for a multi-color flow cytometry experiment.
    """
    # Number of different fluorophores (channels) in the experiment
    num_fluorophores = 5
    
    # 1. Unstained control: for setting baseline fluorescence
    unstained_controls = 1
    
    # 2. Single-stain controls: one for each fluorophore for compensation
    compensation_controls = num_fluorophores
    
    # 3. Fluorescence Minus One (FMO) controls: one for each fluorophore for accurate gating
    fmo_controls = num_fluorophores
    
    # Calculate the total number of control tubes
    total_controls = unstained_controls + compensation_controls + fmo_controls
    
    print("To determine the total number of essential controls, we sum the requirements for each type:")
    print(f"- Unstained Controls: {unstained_controls}")
    print(f"- Single-Stain Compensation Controls: {compensation_controls} (one for each of the {num_fluorophores} fluorophores)")
    print(f"- Fluorescence Minus One (FMO) Controls: {fmo_controls} (one for each of the {num_fluorophores} fluorophores)")
    print("\nCalculating the total:")
    print(f"Total Controls = {unstained_controls} (Unstained) + {compensation_controls} (Compensation) + {fmo_controls} (FMO) = {total_controls}")

calculate_flow_controls()
<<<11>>>