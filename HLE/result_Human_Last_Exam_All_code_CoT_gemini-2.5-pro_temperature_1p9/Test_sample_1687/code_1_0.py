def analyze_blood_loss(initial_hgb, lowest_hgb):
    """
    This function calculates and prints the hemoglobin drop to quantify blood loss.
    The clinical context points to a specific diagnosis based on this data.
    """
    
    # Calculate the drop in hemoglobin
    hemoglobin_drop = initial_hgb - lowest_hgb
    
    # Explain the significance of the calculation
    print("To assess the severity of bleeding, we can calculate the hemoglobin drop.")
    print("The patient's initial hemoglobin was 11.7 g/dL and the lowest recorded level was 6.5 g/dL.")
    
    # Print the equation and result as requested
    print("\nHere is the calculation:")
    print(f"{initial_hgb} g/dL (Initial) - {lowest_hgb} g/dL (Lowest) = {hemoglobin_drop:.1f} g/dL (Total Drop)")
    
    print("\nConclusion based on this data:")
    print("A drop of this magnitude indicates severe internal hemorrhage. This, combined with the location of pain (left upper quadrant and left shoulder), makes splenic laceration the most probable diagnosis following a difficult colonoscopy.")

# Run the analysis with the values from the case study
analyze_blood_loss(initial_hgb=11.7, lowest_hgb=6.5)