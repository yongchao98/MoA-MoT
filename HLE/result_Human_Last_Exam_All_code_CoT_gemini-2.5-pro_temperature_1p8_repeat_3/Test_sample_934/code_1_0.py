import pandas as pd

def diagnose_esophageal_condition():
    """
    Analyzes patient data against known diagnostic criteria for esophageal conditions
    to determine the most likely diagnosis.
    """

    # Patient Profile
    patient_data = {
        'risk_factors': ['heavy_smoking', 'alcohol_use_disorder'],
        'symptoms': ['substernal_chest_pain', 'odynophagia'],
        'labs': ['inflammation'],
        'imaging': ['wall_thickening_narrowing'],
        'endoscopy': ['no_ulcers_plaques_erythema']
    }

    # Diagnostic Criteria
    diagnoses = {
        'A. Streptococcal esophagitis': {
            'strong_risk_factors': ['immunocompromised'],
            'endoscopy_signs': ['exudates', 'ulcers', 'erythema'],
            'common_symptoms': ['fever', 'odynophagia']
        },
        'B. Esophageal adenocarcinoma': {
            'strong_risk_factors': ['chronic_gerd', 'barretts_esophagus', 'obesity'],
            'moderate_risk_factors': ['smoking'],
            'imaging_signs': ['wall_thickening_narrowing', 'mass'],
            'endoscopy_signs': ['mass', 'ulcer', 'nodularity']
        },
        'C. Esophageal squamous cell carcinoma': {
            'strong_risk_factors': ['heavy_smoking', 'alcohol_use_disorder'],
            'imaging_signs': ['wall_thickening_narrowing', 'mass'],
            'symptoms': ['odynophagia', 'chest_pain'],
            # Can present as a submucosal tumor with normal overlying mucosa
            'endoscopy_signs': ['mass', 'ulcer', 'normal_mucosa_possible']
        },
        'D. GERD': {
            'strong_risk_factors': ['obesity', 'hiatal_hernia'],
            'moderate_risk_factors': ['smoking', 'alcohol_use_disorder'],
            'endoscopy_signs': ['erythema', 'erosions', 'strictures']
        },
        'E. Herpes esophagitis': {
            'strong_risk_factors': ['immunocompromised'],
            'endoscopy_signs': ['punched_out_ulcers']
        }
    }

    scores = {name: 0 for name in diagnoses.keys()}
    reasoning = {name: [] for name in diagnoses.keys()}

    print("Evaluating patient findings against potential diagnoses...\n")

    for name, criteria in diagnoses.items():
        # Score risk factors
        if 'strong_risk_factors' in criteria:
            for rf in patient_data['risk_factors']:
                if rf in criteria['strong_risk_factors']:
                    scores[name] += 2
                    reasoning[name].append(f"[+2] Matches strong risk factor: {rf.replace('_', ' ')}")
        if 'moderate_risk_factors' in criteria:
            for rf in patient_data['risk_factors']:
                if rf in criteria['moderate_risk_factors']:
                    scores[name] += 1
                    reasoning[name].append(f"[+1] Matches moderate risk factor: {rf.replace('_', ' ')}")
        
        # Score imaging and labs
        if 'imaging' in patient_data and patient_data['imaging'][0] in criteria.get('imaging_signs', []):
            scores[name] += 2
            reasoning[name].append(f"[+2] Imaging findings (wall thickening) are consistent.")
            
        if 'labs' in patient_data and 'inflammation' in criteria.get('symptoms', []) or name.endswith('carcinoma'):
             if 'inflammation' in patient_data['labs']:
                scores[name] +=1
                reasoning[name].append(f"[+1] Lab findings (inflammation) are consistent.")


        # Score endoscopy findings (crucial step)
        if patient_data['endoscopy'][0] == 'no_ulcers_plaques_erythema':
            # This finding strongly argues against diagnoses defined by these features
            if 'punched_out_ulcers' in criteria.get('endoscopy_signs', []) or \
               'exudates' in criteria.get('endoscopy_signs', []) or \
               'erythema' in criteria.get('endoscopy_signs', []):
                scores[name] -= 3
                reasoning[name].append(f"[-3] Contradictory Endoscopy: Patient lacks expected {criteria['endoscopy_signs'][0].replace('_', ' ')}.")
            # This finding is possible, though not classic, for SCC (submucosal tumor)
            elif 'normal_mucosa_possible' in criteria.get('endoscopy_signs', []):
                scores[name] += 1
                reasoning[name].append(f"[+1] Consistent Endoscopy: A normal mucosal surface is possible.")


    # Determine the most likely diagnosis
    most_likely_diagnosis = max(scores, key=scores.get)

    print("--- Scoring Results ---")
    for name in sorted(scores.keys()):
        print(f"\nDiagnosis: {name}")
        print(f"Final Score: {scores[name]}")
        for reason in reasoning[name]:
            print(f"  - {reason}")
    
    print("\n--- Conclusion ---")
    print(f"The patient has major risk factors (heavy smoking and alcohol use) for Esophageal Squamous Cell Carcinoma.")
    print("Imaging shows wall thickening, consistent with a tumor, while the normal-appearing endoscopy suggests the tumor is submucosal.")
    print("This combination of findings makes other diagnoses less likely.")
    print(f"\nMost Likely Diagnosis: {most_likely_diagnosis}")

diagnose_esophageal_condition()