import sys

def solve_spinal_injury_case():
    """
    Analyzes a clinical vignette to determine the location of a spinal cord injury.
    """
    
    # Step 1: Define the key neurological findings from the vignette.
    patient_findings = {
        "motor_deficit": {"side": "right", "body_part": "leg"},
        "proprioception_vibration_deficit": {"side": "right", "body_part": "leg"},
        "pain_temperature_deficit": {"side": "left", "sensory_level_landmark": "umbilicus"}
    }
    
    # Step 2: Define a map of key dermatome landmarks to their spinal cord levels.
    dermatome_map = {
        "nipple_line": "T4",
        "xiphoid_process": "T6",
        "umbilicus": "T10",
        "inguinal_region": "T12"
    }
    
    # Step 3: Analyze the pattern of deficits to identify the syndrome.
    syndrome_name = "Not identified"
    # Check for Brown-Séquard Syndrome pattern
    if (patient_findings["motor_deficit"]["side"] == patient_findings["proprioception_vibration_deficit"]["side"] and
        patient_findings["motor_deficit"]["side"] != patient_findings["pain_temperature_deficit"]["side"]):
        syndrome_name = "Brown-Séquard Syndrome (Spinal Cord Hemisection)"
        
    # Step 4: Determine the spinal level based on the sensory deficit landmark.
    landmark = patient_findings["pain_temperature_deficit"]["sensory_level_landmark"]
    injury_level = dermatome_map.get(landmark, "Unknown")
    
    # Step 5: Print the detailed reasoning and the final conclusion.
    print("Clinical Reasoning:")
    print("1. The patient presents with a classic triad of symptoms:")
    print(f"   - Ipsilateral (Right) Motor Weakness (Corticospinal Tract lesion)")
    print(f"   - Ipsilateral (Right) Loss of Proprioception/Vibration (Dorsal Column lesion)")
    print(f"   - Contralateral (Left) Loss of Pain/Temperature (Spinothalamic Tract lesion)")
    print("\n2. This constellation of findings is characteristic of Brown-Séquard Syndrome, which is caused by a hemisection of the spinal cord.")
    print("\n3. To determine the level of the injury, we use the sensory level for pain and temperature.")
    print(f"   - The sensory loss begins at the '{landmark.title()}'.")
    print(f"   - The dermatome corresponding to the {landmark} is {injury_level}.")
    
    print("\nConclusion:")
    print(f"The location of the patient's injury is at the {injury_level} spinal level.")
    
    # Step 6: Format the final answer as requested.
    # The answer choice corresponding to T10 is H.
    sys.stdout.write("<<<H>>>\n")

# Execute the function to solve the case.
solve_spinal_injury_case()