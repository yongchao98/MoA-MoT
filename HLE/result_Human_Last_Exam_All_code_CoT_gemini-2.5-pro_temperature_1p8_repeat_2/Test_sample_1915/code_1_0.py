def solve_medical_case():
    """
    This function analyzes the patient's symptoms to determine the location of the anatomical lesion.
    """
    patient_symptoms = {
        "No pupillary light reflex (right eye)": "Points to Cranial Nerve III (Oculomotor) palsy, which carries parasympathetic fibers for pupil constriction.",
        "Inability to adduct (right eye)": "Caused by paralysis of the Medial Rectus muscle, innervated by Cranial Nerve III.",
        "Inability to elevate (right eye)": "Caused by paralysis of the Superior Rectus and Inferior Oblique muscles, both innervated by Cranial Nerve III.",
        "Inability to depress (right eye)": "Caused by paralysis of the Inferior Rectus muscle (innervated by CN III) and Superior Oblique muscle (innervated by CN IV). The collection of symptoms strongly points to a CN III lesion."
    }

    print("Analyzing the Patient's Symptoms:\n")
    for symptom, explanation in patient_symptoms.items():
        print(f"- {symptom}: {explanation}")

    print("\n---")
    print("Conclusion from Symptoms: The combination of these deficits indicates a complete right oculomotor nerve (Cranial Nerve III) palsy.")
    print("---\n")
    
    answer_choices = {
        "A": "Cranial nerve VI: Originates in the Pons. A lesion would cause failure of abduction (outward eye movement).",
        "B": "Cranial nerve VII: Originates in the Pons. A lesion affects muscles of facial expression.",
        "C": "Reticular formation: A diffuse network involved in arousal, not the primary location for these specific motor nuclei.",
        "D": "Medulla oblongata: Contains nuclei for CN IX, X, XI, XII.",
        "E": "Midbrain: Contains the nucleus for Cranial Nerve III (Oculomotor) and Cranial Nerve IV (Trochlear)."
    }

    print("Evaluating Anatomical Locations:\n")
    for choice, description in answer_choices.items():
        print(f"- {choice}. {description}")
    
    print("\n---")
    print("Final Reasoning: A complete Cranial Nerve III palsy is explained by a lesion affecting its nucleus or the nerve fascicles as they exit the brainstem. The nucleus for CN III is located in the Midbrain.")
    print("Therefore, the correct answer is E.")

solve_medical_case()
print("\n<<<E>>>")