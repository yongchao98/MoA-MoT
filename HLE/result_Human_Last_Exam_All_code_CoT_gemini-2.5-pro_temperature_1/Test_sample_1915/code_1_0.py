def solve_clinical_case():
    """
    Analyzes clinical symptoms to identify the damaged anatomical structure.
    """
    # Patient's key symptoms
    symptoms = {
        "Right Eye - Pupillary Light Reflex": "Absent",
        "Right Eye - Adduction (inward move)": "Unable",
        "Right Eye - Depression (downward move)": "Unable",
        "Right Eye - Elevation (upward move)": "Unable"
    }

    # Map of relevant cranial nerve functions
    cn_functions = {
        "CN III (Oculomotor)": ["Pupillary Constriction", "Adduction", "Depression", "Elevation"],
        "CN IV (Trochlear)": ["Depression when adducted"],
        "CN VI (Abducens)": ["Abduction (outward move)"],
        "CN VII (Facial)": ["Facial expression, Eyelid closure"]
    }

    # Map of relevant cranial nerve origins in the brainstem
    cn_origins = {
        "Midbrain": ["CN III", "CN IV"],
        "Pons": ["CN VI", "CN VII"],
        "Medulla Oblongata": ["other CNs"]
    }

    print("Step 1: Identify the primary deficits from the clinical presentation.")
    for symptom, status in symptoms.items():
        print(f"- {symptom}: {status}")

    print("\nStep 2: Determine the affected Cranial Nerve.")
    print("The combination of deficits (pupil reflex, adduction, depression, elevation) overwhelmingly points to a single nerve.")
    affected_cn = "CN III (Oculomotor)"
    print(f"The functions lost correspond directly to the functions of the {affected_cn} nerve.")

    print("\nStep 3: Locate the anatomical origin of the affected nerve.")
    affected_location = ""
    for location, nerves in cn_origins.items():
        if "CN III" in nerves:
            affected_location = location
            break
    
    print(f"The nucleus and initial pathway of the {affected_cn} are located in the {affected_location}.")
    
    print("\nConclusion: The patient's presentation is explained by damage to the Midbrain.")

solve_clinical_case()