def analyze_patient_condition():
    """
    Analyzes the patient's key physical limitation and calculates the deficit.
    """
    # The primary issue is the patient's new inability to ambulate.
    # A key finding is the resistance to knee extension, indicating spasticity.
    # This spasticity is a direct mechanical block to walking.
    # We will quantify this limitation.

    # Data from the clinical vignette
    normal_full_extension_angle = 180  # degrees
    patient_max_extension_angle = 165 # degrees

    # Calculate the knee extension deficit
    extension_deficit = normal_full_extension_angle - patient_max_extension_angle

    # Output the clinical reasoning and the calculation
    print("The patient's key problem is an inability to walk, a decline from his baseline.")
    print("The most specific physical finding explaining this is resistance to left knee extension, a sign of post-stroke spasticity.")
    print("This spasticity is a direct mechanical impediment to walking and effective physical therapy.")
    print("\nWe can quantify this physical limitation:")
    print(f"A normal knee extends to: {normal_full_extension_angle} degrees")
    print(f"This patient's knee extends to: {patient_max_extension_angle} degrees")
    
    print("\nCalculating the extension deficit:")
    print(f"{normal_full_extension_angle} (Normal) - {patient_max_extension_angle} (Patient) = {extension_deficit} (Deficit)")
    
    print(f"\nThis {extension_deficit}-degree extension deficit is the most direct and treatable reason for his inability to ambulate.")
    print("Therefore, the most appropriate next step is to address the spasticity.")

analyze_patient_condition()