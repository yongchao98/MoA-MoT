def generate_denture_plan(patient_info):
    """
    Generates a sample dental treatment plan for a removable partial denture based on patient data.

    This is for informational purposes only and is not medical advice.
    """

    # Extracting data from the input dictionary
    lost_teeth = patient_info["lost_teeth"]
    snb_angle = patient_info["snb_angle"]

    # Define abutment teeth based on the location of lost teeth
    # The teeth adjacent to the gap are chosen as anchors.
    mesial_abutment = "Left Central Incisor"
    distal_abutment = "Left Second Premolar"

    print("--- Prosthetic Treatment Plan: Removable Partial Denture ---")
    print("\nThis plan is formulated after the initial healing phase and medical stabilization.")
    print(f"The patient has lost the following teeth: {', '.join(lost_teeth)}.")

    print("\n1. Type of Denture:")
    print("   - A Cast Partial Denture (CPD) is recommended.")
    print("   - Rationale: It provides superior strength, stability, and support compared to a flexible or all-acrylic denture, which is crucial for replacing multiple teeth.")

    print("\n2. Material Selection:")
    print("   - Framework: Cobalt-Chromium (Co-Cr) alloy. This material is rigid, lightweight, and biocompatible.")
    print("   - Base: High-impact acrylic resin to hold the artificial teeth and sit on the gums.")
    print("   - Artificial Teeth: Acrylic or Composite resin teeth, selected to match the patient's existing teeth.")

    print("\n3. Abutment Teeth and Design Considerations:")
    print(f"   - The primary abutment (support) teeth will be the '{mesial_abutment}' and the '{distal_abutment}'.")
    print("   - Reason for Choice: These teeth are directly adjacent to the edentulous space and can provide the necessary retention and support for the denture via clasps.")
    print("   - Prerequisite: These abutment teeth must be periodontally healthy and free of decay before being used.")
    print(f"\n   - Occlusal Consideration: The skeletal profile (SNB angle of {snb_angle} degrees) must be carefully considered when setting the artificial teeth to ensure a stable and non-traumatic bite, preventing further dental issues.")


# Patient-specific information provided in the prompt
patient_data = {
    "lost_teeth": ["Left Lateral Incisor", "Left Canine", "Left First Premolar"],
    "snb_angle": 88,
    "hbA1c": 7.5,
    "time_since_trauma_hours": 28
}

generate_denture_plan(patient_data)
