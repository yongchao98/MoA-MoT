def explain_diagnosis():
    """
    This function explains the reasoning behind the medical diagnosis.
    """
    patient_profile = {
        "Age": 67,
        "Sex": "Female",
        "Symptoms": ["low-grade fever", "weight loss", "fatigue", "diarrhea (1 month)", "acute RLQ abdominal pain"],
        "Medical History": ["uveitis", "arthritis"],
        "Exam Findings": ["tenderness and guarding in right hemiabdomen"],
        "Lab Findings": {
            "WBC": 13000,
            "Fecal Occult Blood": "Positive"
        },
        "Imaging (CT)": ["marked wall thickening of the ileocecal region", "surrounding mesenteric fat stranding"]
    }

    print("--- Medical Reasoning ---")
    print(f"The patient is a {patient_profile['Age']}-year-old female presenting with a month-long history of constitutional symptoms (fever, weight loss) and gastrointestinal issues (diarrhea), culminating in acute right lower quadrant (RLQ) pain.")
    print("\nKey Diagnostic Clues:")
    print("1. Location of Symptoms/Findings: The pain, tenderness, and CT findings are all localized to the ileocecal region (RLQ), the most common site for Crohn's disease.")
    print("2. Extra-intestinal Manifestations: The patient's history of uveitis and arthritis is highly suggestive of an inflammatory bowel disease, particularly Crohn's disease, for which these are classic associated conditions.")
    print("3. Chronic Inflammatory Picture: The timeline of one month, constitutional symptoms, elevated WBC count ({}), and positive fecal occult blood all point to a chronic inflammatory process, not a simple acute infection.".format(patient_profile['Lab Findings']['WBC']))
    print("4. CT Findings: The CT scan demonstrates significant mural (wall) thickening and inflammation in the ileocecal area, which is a hallmark finding of active Crohn's disease.")

    print("\nEvaluating Other Possibilities:")
    print("- Infections (Yersinia, TB): While possible and can mimic the CT findings, they do not explain the history of uveitis and arthritis as well as Crohn's disease does.")
    print("- Other Diagnoses (Volvulus, Ischemia, Torsion): These are typically acute events with different clinical and imaging presentations and are inconsistent with the patient's one-month history.")
    print("- Malignancy (Lymphoma): This is a consideration, but the entire clinical picture, especially the extra-intestinal manifestations, makes an inflammatory condition like Crohn's disease the most likely diagnosis.")

    print("\n--- Conclusion ---")
    print("The constellation of chronic ileocecal inflammation, constitutional symptoms, and classic extra-intestinal manifestations (uveitis, arthritis) makes Crohn's Disease the most probable diagnosis.")

explain_diagnosis()
print("<<<A>>>")