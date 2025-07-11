def solve_medical_case():
    """
    This script analyzes the patient's symptoms to identify the damaged anatomical structure.
    """

    # Patient's key symptoms for the right eye
    symptoms = {
        "Pupillary Light Reflex": "Absent",
        "Adduction (inward movement)": "Unable",
        "Depression (downward movement)": "Unable",
        "Elevation (upward movement)": "Unable"
    }

    # Mapping of eye functions to Cranial Nerves (CN)
    function_to_nerve = {
        "Pupillary Constriction": "CN III (Oculomotor)",
        "Adduction (Medial Rectus)": "CN III (Oculomotor)",
        "Elevation (Superior Rectus, Inferior Oblique)": "CN III (Oculomotor)",
        "Depression (Inferior Rectus)": "CN III (Oculomotor)",
        "Depression (Superior Oblique)": "CN IV (Trochlear)",
        "Abduction (Lateral Rectus)": "CN VI (Abducens)"
    }

    # Location of cranial nerve nuclei in the brainstem
    nuclei_location = {
        "CN III (Oculomotor)": "Midbrain",
        "CN IV (Trochlear)": "Midbrain",
        "CN V (Trigeminal)": "Pons",
        "CN VI (Abducens)": "Pons",
        "CN VII (Facial)": "Pons",
        "CN VIII (Vestibulocochlear)": "Pons/Medulla",
        "CN IX, X, XI, XII": "Medulla Oblongata"
    }

    print("Step 1: Analyzing the patient's symptoms.")
    print("The patient's right eye cannot adduct, elevate, or depress, and shows no pupillary light reflex.")
    print("-" * 30)

    print("Step 2: Correlating symptoms with cranial nerves.")
    print("The combination of these specific eye movement limitations and the lack of pupillary response points directly to a complete palsy of Cranial Nerve III (the oculomotor nerve).")
    print("-" * 30)

    print("Step 3: Locating the origin of Cranial Nerve III.")
    cn3_location = nuclei_location["CN III (Oculomotor)"]
    print(f"The nucleus for Cranial Nerve III is located in the: {cn3_location}.")
    print("-" * 30)
    
    print("Step 4: Evaluating the answer choices.")
    print("A. Cranial nerve VI lesion would affect outward eye movement (abduction).")
    print("B. Cranial nerve VII lesion would affect facial muscles.")
    print("C. Reticular formation is too general.")
    print("D. Medulla oblongata contains nuclei for lower cranial nerves (IX-XII).")
    print("E. The Midbrain contains the nucleus for Cranial Nerve III.")
    print("-" * 30)

    print("Conclusion: The patient's presentation is best explained by damage to the Midbrain.")

solve_medical_case()