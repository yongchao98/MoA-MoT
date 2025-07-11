import pandas as pd

def diagnose_esophageal_condition():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis.
    """
    patient_findings = {
        'age': 53,
        'symptoms': ['substernal_chest_pain', 'odynophagia'],
        'risk_factors': ['heavy_smoking', 'alcohol_use_disorder'],
        'labs': ['elevated_crp', 'leukocytosis'],
        'imaging': ['esophageal_wall_thickening', 'lumen_narrowing'],
        'endoscopy': ['normal_mucosa']
    }

    diagnoses = {
        'A': {'name': 'Streptococcal esophagitis', 'score': 0, 'reasoning': []},
        'B': {'name': 'Esophageal adenocarcinoma', 'score': 0, 'reasoning': []},
        'C': {'name': 'Esophageal squamous cell carcinoma', 'score': 0, 'reasoning': []},
        'D': {'name': 'GERD', 'score': 0, 'reasoning': []},
        'E': {'name': 'Herpes esophagitis', 'score': 0, 'reasoning': []}
    }

    # --- Scoring Logic ---

    # A. Streptococcal esophagitis
    diagnoses['A']['reasoning'].append("INFO: Very rare diagnosis. Typically presents with fever and exudates on endoscopy.")
    if 'fever' not in patient_findings['symptoms']:
        diagnoses['A']['score'] -= 1
        diagnoses['A']['reasoning'].append("[-1] Patient has no fever.")
    if 'normal_mucosa' in patient_findings['endoscopy']:
        diagnoses['A']['score'] -= 3
        diagnoses['A']['reasoning'].append("[-3] Endoscopy is normal, but strep esophagitis would show plaques/exudates.")
        
    # B. Esophageal adenocarcinoma
    diagnoses['B']['reasoning'].append("INFO: Associated with GERD/Barrett's. Major risk factor is GERD, not alcohol.")
    if 'heavy_smoking' in patient_findings['risk_factors']:
        diagnoses['B']['score'] += 1
        diagnoses['B']['reasoning'].append("[+1] Smoking is a moderate risk factor.")
    if 'normal_mucosa' in patient_findings['endoscopy']:
        diagnoses['B']['score'] -= 5
        diagnoses['B']['reasoning'].append("[-5] MAJOR CONTRADICTION: Adenocarcinoma is a mucosal disease and would almost certainly show a visible mass, ulcer, or nodularity on endoscopy. A normal mucosa is highly unlikely.")

    # C. Esophageal squamous cell carcinoma (SCC)
    diagnoses['C']['reasoning'].append("INFO: Strongly associated with smoking and alcohol.")
    if 'heavy_smoking' in patient_findings['risk_factors']:
        diagnoses['C']['score'] += 3
        diagnoses['C']['reasoning'].append("[+3] Patient has a strong risk factor: heavy smoking.")
    if 'alcohol_use_disorder' in patient_findings['risk_factors']:
        diagnoses['C']['score'] += 3
        diagnoses['C']['reasoning'].append("[+3] Patient has a strong risk factor: alcohol use disorder.")
    if 'esophageal_wall_thickening' in patient_findings['imaging']:
        diagnoses['C']['score'] += 2
        diagnoses['C']['reasoning'].append("[+2] Imaging finding of wall thickening is consistent with a tumor.")
    if 'normal_mucosa' in patient_findings['endoscopy']:
        diagnoses['C']['score'] += 5
        diagnoses['C']['reasoning'].append("[+5] KEY FINDING: SCC can grow submucosally (under the surface), explaining why imaging shows a mass effect while the overlying mucosa on endoscopy appears normal. This reconciles the conflicting test results.")

    # D. GERD
    diagnoses['D']['reasoning'].append("INFO: Common cause of chest pain, but severe symptoms or imaging findings usually accompany endoscopic evidence.")
    if 'normal_mucosa' in patient_findings['endoscopy']:
        diagnoses['D']['score'] -= 4
        diagnoses['D']['reasoning'].append("[-4] Patient's 10/10 pain and odynophagia would likely be caused by severe erosive esophagitis, which is absent on endoscopy.")
    if 'esophageal_wall_thickening' in patient_findings['imaging']:
        diagnoses['D']['score'] -= 2
        diagnoses['D']['reasoning'].append("[-2] Significant wall thickening is not typical for GERD.")

    # E. Herpes esophagitis
    diagnoses['E']['reasoning'].append("INFO: Usually occurs in immunocompromised individuals. Presents with classic 'punched-out' ulcers.")
    diagnoses['E']['score'] -= 2
    diagnoses['E']['reasoning'].append("[-2] Patient is not classically immunocompromised.")
    if 'normal_mucosa' in patient_findings['endoscopy']:
        diagnoses['E']['score'] -= 3
        diagnoses['E']['reasoning'].append("[-3] Endoscopy is normal, but herpes esophagitis would show characteristic vesicles or ulcers.")

    # --- Print Results ---
    print("Evaluating diagnoses based on patient findings:\n")
    
    analysis_data = []
    for key, value in diagnoses.items():
        analysis_data.append({
            'Choice': key,
            'Diagnosis': value['name'],
            'Likelihood Score': value['score'],
            'Reasoning': ' '.join(value['reasoning'])
        })
    
    df = pd.DataFrame(analysis_data)
    pd.set_option('display.max_colwidth', None)
    print(df.to_string(index=False))

    most_likely = max(diagnoses, key=lambda k: diagnoses[k]['score'])
    
    print("\n-------------------------------------------------------------")
    print(f"Conclusion: The diagnosis with the highest likelihood score is '{diagnoses[most_likely]['name']}'.")
    print("This is because the patient's major risk factors (smoking, alcohol) strongly point to it, and it is the only condition listed that classically explains the discrepancy between imaging (showing a mass) and endoscopy (showing normal mucosa) due to its potential for submucosal growth.")
    print("-------------------------------------------------------------")


if __name__ == '__main__':
    diagnose_esophageal_condition()