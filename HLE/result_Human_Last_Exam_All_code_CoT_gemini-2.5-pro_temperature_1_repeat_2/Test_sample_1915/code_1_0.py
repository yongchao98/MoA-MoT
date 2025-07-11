def diagnose_neurological_damage():
    """
    Analyzes clinical symptoms to identify the damaged anatomical structure.
    """
    # Patient's symptoms related to the right eye
    symptoms = {
        "Pupillary Light Reflex": "Absent",
        "Adduction (inward movement)": "Paralyzed",
        "Depression (downward movement)": "Paralyzed",
        "Elevation (upward movement)": "Paralyzed"
    }

    # Functions of relevant anatomical structures
    functions = {
        "Cranial Nerve III (Oculomotor)": ["Pupillary constriction", "Adduction", "Elevation", "Depression (major)"],
        "Cranial Nerve IV (Trochlear)": ["Depression (minor)"],
        "Cranial Nerve VI (Abducens)": ["Abduction (outward movement)"],
        "Midbrain": ["Location of Cranial Nerve III nucleus"],
        "Medulla Oblongata": ["Location of Cranial Nerves IX, X, XI, XII nuclei"]
    }

    print("Patient Presentation Analysis:")
    for symptom, status in symptoms.items():
        print(f"- {symptom}: {status}")

    print("\n--- Diagnostic Logic ---")
    print("1. The combination of symptoms (impaired pupillary reflex, adduction, elevation, and depression) points to a complete palsy of the Oculomotor Nerve (Cranial Nerve III).")
    print(f"2. According to neuroanatomy, the nucleus of Cranial Nerve III is located in the Midbrain.")
    print("3. Therefore, a lesion or traumatic damage to the Midbrain is the most likely cause of the patient's presentation.")
    
    print("\n--- Evaluating Answer Choices ---")
    print("A. Cranial nerve VI damage would only cause inability to move the eye outward.")
    print("B. Cranial nerve VII damage would cause facial paralysis, not these eye movement issues.")
    print("C. Reticular formation damage would affect consciousness more globally.")
    print("D. Medulla oblongata contains nuclei for lower cranial nerves, not CN III.")
    print("E. Midbrain contains the CN III nucleus, making it the correct anatomical location for the damage.")

    print("\nFinal Answer:")
    print("The patient's presentation is explained by damage to the Midbrain.")

if __name__ == "__main__":
    diagnose_neurological_damage()
<<<E>>>