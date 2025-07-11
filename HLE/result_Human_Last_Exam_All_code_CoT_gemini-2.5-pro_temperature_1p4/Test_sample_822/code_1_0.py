import pandas as pd

def diagnose_patient_condition():
    """
    This function analyzes a clinical vignette to determine the most likely diagnosis.
    It simulates a diagnostic process by scoring potential diseases against the patient's symptoms.
    """

    # Define patient symptoms and key findings from the clinical vignette
    patient_profile = {
        "Multi-system inflammation (joints, skin)": True,
        "Multiple pulmonary nodules": True,
        "Constitutional symptoms (fatigue, appetite loss)": True,
        "Neurological symptoms (confusion)": True,
        "History of smoking / Age > 60": True,
        "Immunosuppressed state (from steroids)": True,
        "Fatal opportunistic infection / sepsis": True,
        "Poor response to specific antibiotics (Aminoglycosides)": True,
        "History of asbestos exposure (ship building)": True,
    }

    # Define classic features of differential diagnoses
    disease_features = {
        "Granulomatosis with Polyangiitis (GPA)": {
            "Multi-system inflammation (joints, skin)": True,
            "Multiple pulmonary nodules": True,
            "Constitutional symptoms (fatigue, appetite loss)": True,
            "Neurological symptoms (confusion)": True,
            "History of smoking / Age > 60": False, # Not a direct cause, but common demographic
            "Immunosuppressed state (from steroids)": True, # Disease leads to steroid treatment
            "Fatal opportunistic infection / sepsis": True, # Common complication of disease/treatment
            "Poor response to specific antibiotics (Aminoglycosides)": True, # Infection is often opportunistic
            "History of asbestos exposure (ship building)": False,
        },
        "Lung Cancer with Paraneoplastic Syndrome": {
            "Multi-system inflammation (joints, skin)": True, # Can be caused by paraneoplastic syndrome
            "Multiple pulmonary nodules": True,
            "Constitutional symptoms (fatigue, appetite loss)": True,
            "Neurological symptoms (confusion)": True,
            "History of smoking / Age > 60": True,
            "Immunosuppressed state (from steroids)": True, # Often given as part of treatment
            "Fatal opportunistic infection / sepsis": True, # Possible due to immunosuppression
            "Poor response to specific antibiotics (Aminoglycosides)": False,
            "History of asbestos exposure (ship building)": True,
        },
        "Disseminated Nocardiosis": {
            "Multi-system inflammation (joints, skin)": False, # Infection, not primary inflammation
            "Multiple pulmonary nodules": True,
            "Constitutional symptoms (fatigue, appetite loss)": True,
            "Neurological symptoms (confusion)": True, # Brain abscesses are common
            "History of smoking / Age > 60": False,
            "Immunosuppressed state (from steroids)": True, # Key risk factor
            "Fatal opportunistic infection / sepsis": True,
            "Poor response to specific antibiotics (Aminoglycosides)": True, # Requires specific antibiotics like TMP-SMX
            "History of asbestos exposure (ship building)": False,
        }
    }

    # Analyze and score each potential diagnosis
    analysis_results = []
    for disease, features in disease_features.items():
        score = 0
        explained_symptoms = []
        unexplained_symptoms = []
        for symptom, present in patient_profile.items():
            if present and features.get(symptom, False):
                score += 1
                explained_symptoms.append(symptom)
            elif present:
                unexplained_symptoms.append(symptom)
        
        analysis_results.append({
            "Diagnosis": disease,
            "Score": score,
            "Explained": explained_symptoms,
            "Unexplained": unexplained_symptoms
        })
    
    # Sort results by score to find the best fit
    best_match = sorted(analysis_results, key=lambda x: x['Score'], reverse=True)[0]

    # Print the reasoning
    print("Clinical Reasoning Analysis:")
    print("="*30)
    print("Analyzing the patient's symptoms against potential diagnoses...\n")

    for result in analysis_results:
        print(f"Diagnosis considered: {result['Diagnosis']}")
        print(f"Compatibility Score: {result['Score']} out of {len(patient_profile)}")
        print(f"This diagnosis explains the following findings well:")
        for item in result['Explained']:
            print(f"  - {item}")
        if result['Unexplained']:
            print("This diagnosis does not fully explain:")
            for item in result['Unexplained']:
                print(f"  - {item}")
        print("-" * 20)

    print("\nConclusion:")
    print(f"The disease that provides the most comprehensive explanation for the patient's entire clinical course—from the initial multi-system inflammatory symptoms to the final fatal infection following immunosuppressive therapy—is {best_match['Diagnosis']}.")
    
    return best_match['Diagnosis']

# Execute the diagnosis and print the final answer
final_diagnosis = diagnose_patient_condition()
# The final answer format is handled outside the print statements for clarity.
# In a real application, the return value would be used directly.
# print(f"<<<{final_diagnosis}>>>")
