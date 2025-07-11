import sys
from io import StringIO

def analyze_clinical_case():
    """
    This function analyzes a clinical vignette to determine the most likely diagnosis.
    It models the reasoning process by scoring each option based on provided patient data.
    """
    # Step 1: Define the patient's key clinical features from the vignette.
    patient_features = {
        'symptoms': ['substernal chest pain', 'odynophagia'],
        'risk_factors': ['heavy smoker', 'alcohol use disorder'],
        'findings': ['leukocytosis', 'elevated CRP', 'esophageal wall thickening'],
        'key_endoscopy_finding': 'normal' # No ulcers, plaques, or erythema
    }

    # Step 2: Define the characteristics of each possible diagnosis.
    # A 'contradiction' key is used for findings that strongly argue against a diagnosis.
    diagnoses = {
        'A': {
            'name': 'Streptococcal esophagitis',
            'must_have_endoscopy': ['plaques', 'ulcers', 'erythema'],
            'supports': ['odynophagia', 'leukocytosis']
        },
        'B': {
            'name': 'Esophageal adenocarcinoma',
            'must_have_endoscopy': ['mass', 'ulcer', 'nodule'],
            'supports': ['substernal chest pain'],
            'risk_factors': ['GERD']
        },
        'C': {
            'name': 'Esophageal squamous cell carcinoma',
            # This is the key: Infiltrative SCC can have a normal-appearing mucosa on endoscopy.
            'can_have_normal_endoscopy': True,
            'supports': ['substernal chest pain', 'odynophagia', 'leukocytosis', 'esophageal wall thickening'],
            'risk_factors': ['heavy smoker', 'alcohol use disorder']
        },
        'D': {
            'name': 'GERD',
            # Severe GERD causing this pain would show signs of esophagitis.
            'must_have_endoscopy': ['erythema', 'erosions', 'ulcers'],
            'supports': ['heartburn', 'regurgitation'] # These are absent in the patient.
        },
        'E': {
            'name': 'Herpes esophagitis',
            'must_have_endoscopy': ['punched-out ulcers'],
            'supports': ['odynophagia']
        }
    }

    # Step 3: Score each diagnosis based on the patient's features.
    scores = {}
    analysis_log = ""
    # Use a string buffer to capture print output for analysis log
    old_stdout = sys.stdout
    sys.stdout = analysis_buffer = StringIO()
    
    print("--- Diagnostic Reasoning ---")
    
    for key, data in diagnoses.items():
        score = 0
        reasoning = f"\nAnalyzing Diagnosis {key}: {data['name']}"
        
        # Check for strong contradictions
        if 'must_have_endoscopy' in data and patient_features['key_endoscopy_finding'] == 'normal':
            score -= 10
            reasoning += f"\n  - Major Contradiction: This diagnosis typically shows {data['must_have_endoscopy']} on endoscopy, but the patient's was normal. Score -10."

        # Check for strong positive correlation with the key finding
        if data.get('can_have_normal_endoscopy') and patient_features['key_endoscopy_finding'] == 'normal':
            score += 5
            reasoning += "\n  - Key Feature Match: This diagnosis can present with a normal endoscopy (infiltrative type), which matches the patient's case. Score +5."

        # Add points for matching risk factors
        for risk in patient_features['risk_factors']:
            if risk in data.get('risk_factors', []):
                score += 3
                reasoning += f"\n  - Risk Factor Match: Patient has '{risk}'. Score +3."
        
        # Add points for matching symptoms and findings
        for feature_list in ['symptoms', 'findings']:
            for feature in patient_features[feature_list]:
                if feature in data.get('supports', []):
                    score += 1
                    reasoning += f"\n  - Finding Match: Patient has '{feature}'. Score +1."
        
        scores[key] = score
        print(reasoning)
        print(f"  --> Final Score for {key}: {score}")

    # Determine the best diagnosis
    most_likely_diagnosis_key = max(scores, key=scores.get)
    
    # Restore stdout and get the log
    sys.stdout = old_stdout
    analysis_log = analysis_buffer.getvalue()

    # Final Output
    print(analysis_log)
    print("\n--- Conclusion ---")
    print(f"The most likely diagnosis is '{diagnoses[most_likely_diagnosis_key]['name']}' because it aligns best with the patient's major risk factors (smoking, alcohol) and uniquely explains the combination of wall thickening on imaging with a normal endoscopic appearance.")
    
    # Final Answer format
    print("\n<<<C>>>")

# Execute the analysis
analyze_clinical_case()