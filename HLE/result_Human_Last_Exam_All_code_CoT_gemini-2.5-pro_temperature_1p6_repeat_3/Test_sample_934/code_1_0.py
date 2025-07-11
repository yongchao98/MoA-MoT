def diagnose_esophageal_condition():
    """
    This function analyzes a patient's clinical data to determine the most
    likely esophageal diagnosis by assigning scores based on key findings.
    """
    patient_findings = {
        'heavy_smoking': 1,
        'alcohol_use_disorder': 1,
        'odynophagia': 1,
        'imaging_wall_thickening': 1,
        'endoscopy_normal_mucosa': 1,
        'endoscopy_ulcers_or_plaques': 0,
        'history_of_gerd': 0
    }

    diagnoses = {
        'A. Streptococcal esophagitis': {
            'heavy_smoking': 0, 'alcohol_use_disorder': 0, 'odynophagia': 1,
            'imaging_wall_thickening': 0, 'endoscopy_normal_mucosa': 0,
            'endoscopy_ulcers_or_plaques': 1, 'history_of_gerd': 0
        },
        'B. Esophageal adenocarcinoma': {
            'heavy_smoking': 0.5, 'alcohol_use_disorder': 0.5, 'odynophagia': 1,
            'imaging_wall_thickening': 1, 'endoscopy_normal_mucosa': 0,
            'endoscopy_ulcers_or_plaques': 1, 'history_of_gerd': 1
        },
        'C. Esophageal squamous cell carcinoma': {
            'heavy_smoking': 1, 'alcohol_use_disorder': 1, 'odynophagia': 1,
            'imaging_wall_thickening': 1, 'endoscopy_normal_mucosa': 1, # Classic for infiltrative type
            'endoscopy_ulcers_or_plaques': 0, 'history_of_gerd': 0
        },
        'D. GERD': {
            'heavy_smoking': 0.5, 'alcohol_use_disorder': 0.5, 'odynophagia': 0.5,
            'imaging_wall_thickening': 0, 'endoscopy_normal_mucosa': 0,
            'endoscopy_ulcers_or_plaques': 1, 'history_of_gerd': 1
        },
        'E. Herpes esophagitis': {
            'heavy_smoking': 0, 'alcohol_use_disorder': 0, 'odynophagia': 1,
            'imaging_wall_thickening': 0, 'endoscopy_normal_mucosa': 0,
            'endoscopy_ulcers_or_plaques': 1, 'history_of_gerd': 0
        }
    }

    results = {}
    print("Evaluating patient findings against potential diagnoses...\n")

    for diagnosis, criteria in diagnoses.items():
        score = 0
        equation_parts = []
        for finding, patient_value in patient_findings.items():
            # Only score points if the finding is present in the patient
            if patient_value > 0:
                # Add score if the diagnosis matches the patient's finding
                match_score = criteria.get(finding, 0)
                if match_score > 0:
                    score += match_score
                    equation_parts.append(f"({finding}) {match_score}")
        
        results[diagnosis] = score
        print(f"Diagnosis: {diagnosis}")
        equation_str = " + ".join(equation_parts)
        print(f"Score Equation: {equation_str} = {score:.1f}")
        print("-" * 30)

    most_likely_diagnosis = max(results, key=results.get)
    print(f"\nConclusion: The highest score belongs to '{most_likely_diagnosis}', making it the most likely diagnosis.")

if __name__ == '__main__':
    diagnose_esophageal_condition()