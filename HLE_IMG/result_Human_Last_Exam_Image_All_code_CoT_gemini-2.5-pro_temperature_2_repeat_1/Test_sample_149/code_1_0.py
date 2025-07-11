def provide_diagnosis():
    """
    This function analyzes the provided clinical and histological data to provide a diagnosis.
    """
    patient_age = 28
    patient_gender = "female"
    lesion_location = "sella turcica"
    
    # Histopathological findings from the PAS-stained image
    histology = {
        "cell_infiltrate": "Prominent inflammatory infiltrate (lymphocytes)",
        "key_finding": "Presence of fungal organisms",
        "organism_morphology": "Branching, septate hyphae"
    }
    
    # Diagnosis based on synthesis of findings
    # The presence of septate fungal hyphae in a sellar mass is characteristic of an Aspergillus infection,
    # often extending from the adjacent paranasal sinuses.
    diagnosis = "Aspergillosis"
    
    print("Clinical and Histopathological Analysis:")
    print(f"Patient: {patient_age}-year-old {patient_gender}")
    print(f"Location of Mass: {lesion_location}")
    print("Histology:")
    print(f"  - Infiltrate: {histology['cell_infiltrate']}")
    print(f"  - Key Finding: {histology['key_finding']}")
    print(f"  - Organism Morphology: {histology['organism_morphology']}")
    print("-" * 30)
    print(f"Final Diagnosis: {diagnosis}")

provide_diagnosis()