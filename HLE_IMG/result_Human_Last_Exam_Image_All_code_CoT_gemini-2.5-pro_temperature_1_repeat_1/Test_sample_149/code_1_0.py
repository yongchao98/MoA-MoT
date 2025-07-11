def provide_diagnosis():
    """
    Analyzes the clinical and histopathological data to provide a diagnosis.
    """
    
    # Clinical Information
    patient_age = 28
    patient_gender = "female"
    lesion_location = "mass in the sella turcica"

    # Histopathological Findings
    stain = "Periodic acid-Schiff (PAS) stain, highlighting fungal elements."
    hyphae_morphology = "Broad, ribbon-like, pauci-septate (largely non-septate) hyphae."
    branching_pattern = "Irregular, wide-angle branching, often at 90 degrees."
    background_tissue = "Necrosis and a dense inflammatory infiltrate."

    # Synthesis and Conclusion
    conclusion = """The combination of broad, pauci-septate, ribbon-like hyphae with right-angle branching is classic for Mucormycosis.
The location in the sella turcica suggests an invasive fungal process, likely extending from the adjacent sphenoid sinus."""

    # Printing the diagnostic summary
    print("Diagnosis Derivation:")
    print("-" * 20)
    print(f"Patient: {patient_age}-year-old {patient_gender} with a {lesion_location}.")
    print("\nHistopathology shows:")
    print(f"- Stain: {stain}")
    print(f"- Fungal Morphology: {hyphae_morphology}")
    print(f"- Branching Pattern: {branching_pattern}")
    print(f"- Background: {background_tissue}")
    print("\nConclusion:")
    print(conclusion)
    print("\nFinal Diagnosis: Mucormycosis")

# Execute the function to print the diagnosis
provide_diagnosis()