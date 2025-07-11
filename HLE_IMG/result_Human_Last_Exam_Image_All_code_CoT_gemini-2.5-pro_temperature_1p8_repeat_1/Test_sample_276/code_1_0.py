def diagnose_patient():
    """
    Analyzes the patient's clinical data to determine the most likely diagnosis.
    """
    # Patient Data
    patient = {
        "age": 67,
        "sex": "female",
        "symptoms": ["low-grade fever", "weight loss", "fatigue", "diarrhea (1 month)", "acute right abdominal pain"],
        "history": ["gallstones", "uveitis", "arthritis"],
        "exam_findings": ["right hemiabdomen tenderness", "guarding"],
        "labs": {
            "WBC_count": 13000,  # cells/mcL
            "FOBT": "positive"
        },
        "imaging": ["CT scan showing marked circumferential wall thickening of the ileocecal region with fat stranding"]
    }

    # Analysis
    print("Patient Clinical Profile Analysis:")
    print(f"- A {patient['age']}-year-old female presents with a 1-month history of constitutional symptoms (fever, weight loss) and diarrhea.")
    print(f"- She has a significant medical history of uveitis and arthritis.")
    print(f"- Physical exam and imaging localize the pathology to the ileocecal region (right lower quadrant).")
    print(f"- Labs show inflammation (WBC: {patient['labs']['WBC_count']}) and GI bleeding (FOBT: {patient['labs']['FOBT']}).")
    print("- The CT scan confirms active inflammation with bowel wall thickening in the classic location for Crohn's Disease.")

    # Conclusion
    conclusion = "Crohn's Disease (A) is the most likely diagnosis as it cohesively explains the chronic GI symptoms, the constitutional symptoms, the acute presentation, the specific location of inflammation, AND the patient's history of extraintestinal manifestations (uveitis and arthritis)."
    
    print("\nMost Likely Diagnosis:")
    print(conclusion)

# Execute the function to display the analysis.
diagnose_patient()