def calculate_flow_controls():
    """
    Calculates the number of essential controls for a multi-color flow cytometry
    sorting experiment.
    """
    # 1. Define experiment parameters based on the user's question
    fluorophores = ["AF350", "GFP", "PE", "AF647", "AF750"]
    num_channels = len(fluorophores)
    uses_streptavidin_reagent = True  # Based on the mention of "Streptavidin beads"

    # 2. Calculate the number of each essential control type
    
    # Unstained control: For setting baseline autofluorescence
    unstained_control_count = 1
    
    # Single-stain controls: For compensation (one for each fluorophore)
    compensation_control_count = num_channels
    
    # FMO controls: For accurate gating in sorting (one for each fluorophore)
    fmo_control_count = num_channels
    
    # Streptavidin reagent control: To check for non-specific binding of the streptavidin
    streptavidin_reagent_control_count = 1 if uses_streptavidin_reagent else 0
    
    # 3. Sum the counts to get the total
    total_controls = (unstained_control_count + 
                      compensation_control_count + 
                      fmo_control_count + 
                      streptavidin_reagent_control_count)

    # 4. Print the breakdown and the final equation
    print("To perform a robust 5-color sorting experiment with a streptavidin system, you should prepare the following essential controls:")
    print("-" * 40)
    print(f"Unstained Controls: {unstained_control_count}")
    print(f"Single-Stain Compensation Controls: {compensation_control_count} (one for each fluorophore)")
    print(f"Fluorescence Minus One (FMO) Controls: {fmo_control_count} (for accurate gating)")
    print(f"Streptavidin Reagent Control: {streptavidin_reagent_control_count} (for non-specific binding)")
    print("-" * 40)
    
    # Final equation as requested
    print("Total Essential Controls Calculation:")
    print(f"Unstained ({unstained_control_count}) + Compensation ({compensation_control_count}) + FMO ({fmo_control_count}) + Streptavidin ({streptavidin_reagent_control_count}) = {total_controls}")

calculate_flow_controls()
<<<12>>>