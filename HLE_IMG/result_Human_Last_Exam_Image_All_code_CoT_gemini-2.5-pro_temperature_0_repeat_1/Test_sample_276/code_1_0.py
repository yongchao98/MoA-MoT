def solve_medical_case():
    """
    This function analyzes the provided clinical case and determines the most likely diagnosis.
    """
    patient_profile = {
        "age": 67,
        "sex": "female",
        "symptoms": ["low-grade fever", "weight loss", "fatigue", "diarrhea", "acute abdominal pain"],
        "symptom_duration": "1 month",
        "past_medical_history": ["gallstones", "uveitis", "arthritis"],
        "exam_findings": ["right hemiabdomen tenderness", "guarding"],
        "lab_results": {
            "WBC": 13000, # leukocytosis
            "fecal_occult_blood": "positive"
        },
        "imaging": "CT scan showing marked ileocecal wall thickening and surrounding inflammation"
    }

    # Analysis: The key features are the location of inflammation (ileocecal),
    # the chronic nature of the symptoms, and the extraintestinal manifestations
    # (uveitis, arthritis). This constellation of findings is classic for Crohn's Disease.
    # While other conditions like TB or lymphoma can mimic the imaging, the full
    # clinical picture strongly supports Crohn's.

    answer_choices = {
        "A": "Crohn's Disease",
        "B": "Yersinia Colitis",
        "C": "Ileocecal Tuberculosis",
        "D": "Salmonella Enteritis",
        "E": "C. difficile Colitis",
        "F": "Cecal Volvulus",
        "G": "Ischemic Colitis",
        "H": "Pneumoperitoneum",
        "I": "Ovarian Torsion",
        "J": "Celiac Disease",
        "K": "Gastrointestinal Lymphoma"
    }

    most_likely_diagnosis_code = "A"
    
    print(f"The patient's presentation with chronic diarrhea, weight loss, fever, and acute right lower quadrant pain is concerning.")
    print(f"Laboratory findings show leukocytosis (WBC: {patient_profile['lab_results']['WBC']}) and gastrointestinal bleeding (FOBT: {patient_profile['lab_results']['fecal_occult_blood']}).")
    print(f"The CT scan reveals significant inflammation and thickening in the ileocecal region.")
    print(f"Crucially, the patient's history of {patient_profile['past_medical_history'][1]} and {patient_profile['past_medical_history'][2]} are known extraintestinal manifestations of Inflammatory Bowel Disease.")
    print(f"Given that Crohn's Disease classically affects the ileocecal region and is associated with these systemic inflammatory conditions, it is the most likely diagnosis.")
    print(f"The final answer is {most_likely_diagnosis_code}: {answer_choices[most_likely_diagnosis_code]}")

solve_medical_case()