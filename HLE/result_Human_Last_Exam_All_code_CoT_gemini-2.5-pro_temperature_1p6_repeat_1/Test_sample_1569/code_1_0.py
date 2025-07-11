import json

def solve_clinical_case():
    """
    Analyzes a clinical case by scoring potential diagnoses based on
    geography, symptoms, and lab findings.
    """

    # --- 1. Define Disease Profiles ---
    disease_profiles = {
        'A': {
            "name": "Babesia microti",
            "geography": ["Northeast US", "Upper Midwest US"],
            "symptoms": ["Fever", "Headaches", "Myalgia"],
            "lab_notes": "Cross-reactivity with Lyme is possible but less common."
        },
        'B': {
            "name": "Plasmodium",
            "geography": ["Africa", "Asia", "South America"],
            "symptoms": ["Fever", "Headaches", "Myalgia", "Disorientation"],
            "lab_notes": "Not known for Lyme serology cross-reactivity."
        },
        'C': {
            "name": "Borrelia burgdorferi",
            "geography": ["Northeast US", "Upper Midwest US", "Mid-Atlantic US"],
            "symptoms": ["Fever", "Headaches", "Myalgia", "Disorientation", "Heart murmur"],
            "lab_notes": "Negative IgG makes active infection with these symptoms less likely."
        },
        'D': {
            "name": "Ehrlichia",
            "geography": ["Southeastern US", "South-central US", "Oklahoma"],
            "symptoms": ["Fever", "Headaches", "Myalgia", "Disorientation"],
            "lab_notes": "Known to cause false-positive Lyme IgM serology."
        },
        'E': {
            "name": "Rickettsia rickettsii",
            "geography": ["Southeastern US", "Oklahoma"],
            "symptoms": ["Fever", "Headaches", "Myalgia", "Disorientation"],
            "lab_notes": "Serologic cross-reactivity with other pathogens can occur."
        }
    }

    # --- 2. Define Patient Case ---
    patient_case = {
        "geography": "Oklahoma",
        "symptoms": ["Fever", "Headaches", "Myalgia", "Disorientation", "Heart murmur"],
        "labs": "Positive IgM, Negative IgG Lyme"
    }

    # --- 3. Score each diagnosis ---
    print("Evaluating potential diagnoses based on the patient's case...\n")
    scores = {}

    for key, disease in disease_profiles.items():
        score = 0
        score_breakdown = {}

        # Score Geography (Weight: 3)
        geo_score = 3 if patient_case["geography"] in disease["geography"] else -5
        score += geo_score
        score_breakdown["Geography"] = geo_score

        # Score Symptoms (Weight: 1 per match)
        symptom_score = sum(1 for symptom in patient_case["symptoms"] if symptom in disease["symptoms"])
        score += symptom_score
        score_breakdown["Symptoms"] = symptom_score

        # Score Lab context (Weight: 2)
        lab_score = 0
        if "false-positive Lyme" in disease["lab_notes"]:
            lab_score = 2
        elif key == 'C' and "Negative IgG" in patient_case["labs"]:
            # Penalize Lyme because negative IgG argues against it
            lab_score = -2
        score += lab_score
        score_breakdown["Labs"] = lab_score

        scores[key] = {
            "name": disease["name"],
            "total_score": score,
            "breakdown": score_breakdown
        }

    # --- 4. Display Results and reasoning ---
    sorted_scores = sorted(scores.items(), key=lambda item: item[1]['total_score'], reverse=True)

    for key, result in sorted_scores:
        name = result['name']
        total = result['total_score']
        breakdown = result['breakdown']
        
        # This part fulfills the prompt's requirement to output each number in the final equation
        geo_num = breakdown['Geography']
        sxs_num = breakdown['Symptoms']
        lab_num = breakdown['Labs']
        
        print(f"Diagnosis: ({key}) {name}")
        print(f"Final Score: {total}")
        print(f"Score Equation: (Geography) {geo_num} + (Symptoms) {sxs_num} + (Labs) {lab_num} = {total}")
        print("-" * 30)

    best_fit_key = sorted_scores[0][0]
    print(f"\nConclusion: The highest score belongs to ({best_fit_key}) {scores[best_fit_key]['name']}, making it the most likely diagnosis.")

if __name__ == '__main__':
    solve_clinical_case()
