def diagnose_neurological_lesion():
    """
    Analyzes patient symptoms to identify the location of a neurological lesion.
    """
    symptoms = {
        "Pupillary Light Reflex": "Absent",
        "Adduction (inward movement)": "Unable",
        "Depression (downward movement)": "Unable",
        "Elevation (upward movement)": "Unable"
    }

    cranial_nerve_functions = {
        "CN III (Oculomotor)": ["Pupillary Constriction", "Adduction", "Elevation", "Depression"],
        "CN IV (Trochlear)": ["Depression (via Superior Oblique)"],
        "CN VI (Abducens)": ["Abduction (outward movement)"]
    }

    brainstem_nuclei_location = {
        "Midbrain": ["CN III", "CN IV"],
        "Pons": ["CN V", "CN VI", "CN VII", "CN VIII"],
        "Medulla Oblongata": ["CN IX", "CN X", "CN XI", "CN XII"]
    }

    print("Step 1: Analyzing patient's right eye deficits.")
    affected_functions = []
    for symptom, status in symptoms.items():
        print(f"- {symptom}: {status}")
        if "Pupil" in symptom:
            affected_functions.append("Pupillary Constriction")
        else:
            affected_functions.append(symptom.split(' ')[0])

    print("\nStep 2: Identifying the affected cranial nerve.")
    print("The combination of deficits in pupillary constriction, adduction, elevation, and depression points to a complete palsy of Cranial Nerve III (Oculomotor Nerve).")

    print("\nStep 3: Locating the lesion based on neuroanatomy.")
    location = ""
    for region, nerves in brainstem_nuclei_location.items():
        if "CN III" in nerves:
            location = region
            break
    
    print(f"The nucleus of Cranial Nerve III is located in the {location}.")
    print("Therefore, damage from the patient's trauma and stroke most likely occurred in this area.")

    print("\nConclusion: The patient's presentation is explained by damage to the Midbrain.")

diagnose_neurological_lesion()