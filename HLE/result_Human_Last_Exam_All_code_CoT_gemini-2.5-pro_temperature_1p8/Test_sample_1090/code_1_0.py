def analyze_patient_mobility_issue():
    """
    This script analyzes the key physical finding preventing the patient's ambulation.
    """
    
    # Define the parameters for knee extension based on the clinical case.
    # Full extension of the knee joint is 180 degrees.
    full_knee_extension = 180
    
    # The patient's paretic left knee cannot be extended beyond 165 degrees due to resistance.
    patients_knee_extension = 165
    
    # Calculate the patient's knee extension deficit.
    extension_deficit = full_knee_extension - patients_knee_extension
    
    print("Analyzing the patient's inability to ambulate:")
    print("A key finding from the physical exam is a limitation in the left leg, which was affected by a previous stroke.")
    print("\nTo quantify this limitation, we calculate the knee extension deficit:")
    # We output each number in the final equation as requested.
    print(f"Normal Full Extension ({full_knee_extension}°) - Patient's Maximum Extension ({patients_knee_extension}°) = Deficit ({extension_deficit}°)")
    
    print(f"\nConclusion:")
    print(f"A {extension_deficit}-degree extension deficit is significant.")
    print("This indicates severe spasticity or a developing contracture, which acts as a direct mechanical block to standing and walking.")
    print("While the patient is deconditioned, this specific, treatable issue is the most critical barrier to address to restore his mobility.")

analyze_patient_mobility_issue()