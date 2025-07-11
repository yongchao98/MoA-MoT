def diagnosis_reasoning():
    """
    This function provides a step-by-step analysis of the clinical case to reach a diagnosis.
    """
    patient_profile = {
        "Age": 67,
        "Sex": "Female",
        "Symptoms": ["low-grade fever", "weight loss", "fatigue", "diarrhea (1 month)", "acute abdominal pain"],
        "History": ["gallstones", "uveitis", "arthritis"],
        "Exam Findings": ["right hemiabdomen tenderness", "guarding"],
        "Lab Results": {
            "WBC Count": 13000,  # cells/mcL
            "Fecal Occult Blood": "Positive"
        },
        "Imaging": "CT scan showing marked circumferential wall thickening of the ileocecal region with surrounding inflammation."
    }

    print("Clinical Case Analysis:")
    print("-" * 25)
    print(f"Patient: {patient_profile['Age']}-year-old {patient_profile['Sex']}")
    print(f"Presentation: A 1-month history of constitutional symptoms ({', '.join(patient_profile['Symptoms'][0:4])}) and new acute pain.")
    print(f"Key History: The patient has a history of {patient_profile['History'][1]} and {patient_profile['History'][2]}, which are classic extraintestinal manifestations of Inflammatory Bowel Disease (IBD).")
    print(f"Exam & Labs: Right lower quadrant findings, elevated WBC ({patient_profile['Lab Results']['WBC Count']}), and positive fecal occult blood all point to active inflammation and bleeding in the GI tract.")
    print(f"Imaging Findings: The CT scan localizes the problem to the ileocecal region, showing severe inflammation. This is the most common site for Crohn's Disease.")
    print("\nDifferential Diagnosis Evaluation:")
    print(" - Crohn's Disease (A): Excellent fit. Matches the location (ileocecal), chronicity, extraintestinal manifestations, and CT findings.")
    print(" - Ileocecal TB (C): A strong mimic, but the history of uveitis and arthritis makes Crohn's more likely.")
    print(" - Yersinia (B), Salmonella (D), C. diff (E): Less likely due to the chronic nature and specific extraintestinal symptoms.")
    print(" - Other options (F, G, H, I, J, K): Do not fit the combination of chronic and acute symptoms, extraintestinal manifestations, and specific CT findings.")
    print("\nConclusion:")
    print("The combination of chronic ileocecal inflammation, constitutional symptoms, and a history of uveitis and arthritis makes Crohn's Disease the most likely diagnosis.")

diagnosis_reasoning()