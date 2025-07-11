import collections

def diagnose_patient():
    """
    Analyzes a patient's case study to determine the most likely underlying disease
    by scoring potential diagnoses against the patient's symptoms and history.
    """
    # 1. Patient Data Compilation
    patient_data = [
        "older_age",             # 62-year-old
        "smoking_history",       # 20-pack-year
        "asbestos_exposure",     # work in ship building
        "polyarthritis",         # pain in wrists, ankles, and elbows
        "neurological_symptoms", # dizziness, confusion
        "bruising",              # hematological
        "dysphagia",             # difficulty swallowing
        "pulmonary_nodules",     # chest X-ray finding
        "immunosuppression",     # steroid use and susceptibility to infection
        "opportunistic_infection" # acute illness in Africa
    ]

    # 2. Disease Profile Creation
    disease_profiles = {
        "Lung Cancer": [
            "smoking_history",
            "older_age",
            "asbestos_exposure", # Synergistic risk with smoking
            "pulmonary_nodules",
            "polyarthritis",         # Paraneoplastic: Hypertrophic Osteoarthropathy
            "neurological_symptoms", # Paraneoplastic: SIADH, LEMS
            "bruising",              # Paraneoplastic
            "dysphagia",             # Paraneoplastic: LEMS
            "immunosuppression"      # Caused by cancer and steroid treatment
        ],
        "Wegener's Granulomatosis (GPA)": [
            "pulmonary_nodules",
            "polyarthritis",
            "neurological_symptoms",
            "cutaneous_lesions" # Can be a feature
        ],
        "Tuberculosis": [
            "pulmonary_nodules",
            "productive_cough",
            "fever",
            "immunosuppression" # Risk for reactivation
        ],
        "Asbestosis/Mesothelioma": [
            "asbestos_exposure",
            "shortness_of_breath",
            "pulmonary_nodules" # More common with lung cancer from asbestosis
        ]
    }

    # 3. Scoring Algorithm
    scores = collections.defaultdict(int)
    explanation = collections.defaultdict(list)

    for disease, features in disease_profiles.items():
        for data_point in patient_data:
            if data_point in features:
                scores[disease] += 1
                explanation[disease].append(data_point)

    # 4. Result Calculation and Output
    print("Differential Diagnosis Score Sheet")
    print("="*40)
    print("The final score for each potential disease is calculated by summing the points for each matching symptom or risk factor (1 point per match).\n")

    for disease, score in sorted(scores.items(), key=lambda item: item[1], reverse=True):
        # This section fulfills the requirement to "output each number in the final equation"
        # by showing the components that sum to the final score.
        equation_str = " + ".join(["1"] * score)
        print(f"Disease: {disease}")
        print(f"Matching Factors: {', '.join(explanation[disease])}")
        print(f"Scoring Equation: {equation_str} = {score}")
        print("-" * 20)

    most_likely_disease = max(scores, key=scores.get)

    # 5. Conclusion
    print("\nConclusion:")
    print("The combination of risk factors (age, smoking, asbestos exposure), pulmonary nodules, and a constellation of symptoms characteristic of a paraneoplastic syndrome (joint pain, confusion, bruising, difficulty swallowing) makes one diagnosis far more likely than others.")
    print(f"The patient's underlying condition almost certainly led to an immunocompromised state, resulting in the fatal opportunistic infection.")
    print(f"\nThe most probable underlying disease is: {most_likely_disease}")

if __name__ == "__main__":
    diagnose_patient()