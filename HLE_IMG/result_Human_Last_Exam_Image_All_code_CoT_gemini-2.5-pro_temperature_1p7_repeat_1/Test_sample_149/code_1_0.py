def provide_diagnosis():
    """
    Analyzes the clinical and histological findings to provide a diagnosis.
    """
    # Clinical Information
    patient_age = 28
    patient_sex = "female"
    location = "mass in the sella turcica"

    # Histological Findings
    stain = "Periodic acid-Schiff (PAS)"
    fungal_morphology = [
        "Broad, ribbon-like hyphae",
        "Pauciseptate (sparsely septate)",
        "Right-angle or obtuse-angle branching",
        "Invasion of tissue with associated inflammation"
    ]

    # Diagnosis
    diagnosis = "Mucormycosis"
    
    # Explanation
    print(f"Patient Information: A {patient_age}-year-old female presents with a {location}.")
    print("\nHistological Analysis:")
    print(f"The PAS stain reveals fungal organisms with the following key features:")
    for feature in fungal_morphology:
        print(f"- {feature}")
    
    print("\nConclusion:")
    print("These morphological features—specifically the broad, pauciseptate hyphae with right-angle branching—are characteristic of fungi in the order Mucorales.")
    print(f"Therefore, the final diagnosis is: {diagnosis}")

if __name__ == "__main__":
    provide_diagnosis()