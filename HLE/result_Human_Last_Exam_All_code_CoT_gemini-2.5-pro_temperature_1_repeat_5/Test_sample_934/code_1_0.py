import pandas as pd

def diagnose_esophageal_condition():
    """
    This function analyzes clinical findings to determine the most likely esophageal diagnosis.
    It uses a scoring system where positive points are awarded for findings that support a diagnosis,
    and negative points for findings that contradict it.
    """
    
    # Clinical findings from the case vignette
    findings = {
        'Major Risk Factors (Smoking/Alcohol)': 2, # Strong presence
        'Symptoms (Severe Odynophagia/Chest Pain)': 1, # Present
        'Imaging (Wall Thickening/Narrowing)': 2, # Present and significant
        'Normal Endoscopy (No Mucosal Lesions)': 2, # A key, differentiating finding
        'Inflammatory Labs (Leukocytosis/CRP)': 1 # Present
    }

    # Scoring matrix for each diagnosis based on clinical knowledge
    # Score = How much a finding supports (+) or refutes (-) a diagnosis
    scoring_matrix = {
        'Streptococcal esophagitis': {
            'Major Risk Factors (Smoking/Alcohol)': 0,
            'Symptoms (Severe Odynophagia/Chest Pain)': 1,
            'Imaging (Wall Thickening/Narrowing)': 1,
            'Normal Endoscopy (No Mucosal Lesions)': -2, # Expects exudates/ulcers
            'Inflammatory Labs (Leukocytosis/CRP)': 1
        },
        'Esophageal adenocarcinoma': {
            'Major Risk Factors (Smoking/Alcohol)': 1, # Smoking is a risk factor, but GERD is primary
            'Symptoms (Severe Odynophagia/Chest Pain)': 1,
            'Imaging (Wall Thickening/Narrowing)': 2,
            'Normal Endoscopy (No Mucosal Lesions)': -2, # Expects visible mass/ulcer/Barrett's
            'Inflammatory Labs (Leukocytosis/CRP)': 1
        },
        'Esophageal squamous cell carcinoma': {
            'Major Risk Factors (Smoking/Alcohol)': 2, # Classic, strong risk factors
            'Symptoms (Severe Odynophagia/Chest Pain)': 1,
            'Imaging (Wall Thickening/Narrowing)': 2,
            'Normal Endoscopy (No Mucosal Lesions)': 2, # Consistent with infiltrative/submucosal cancer
            'Inflammatory Labs (Leukocytosis/CRP)': 1 # Consistent with paraneoplastic syndrome
        },
        'GERD': {
            'Major Risk Factors (Smoking/Alcohol)': 1, # Can worsen GERD
            'Symptoms (Severe Odynophagia/Chest Pain)': 0, # Chest pain fits, but severe odynophagia is less typical without erosions
            'Imaging (Wall Thickening/Narrowing)': -1, # Not a typical finding for GERD
            'Normal Endoscopy (No Mucosal Lesions)': -2, # Severe symptoms should cause visible esophagitis
            'Inflammatory Labs (Leukocytosis/CRP)': -1 # Not typical for uncomplicated GERD
        },
        'Herpes esophagitis': {
            'Major Risk Factors (Smoking/Alcohol)': 0,
            'Symptoms (Severe Odynophagia/Chest Pain)': 1,
            'Imaging (Wall Thickening/Narrowing)': 0,
            'Normal Endoscopy (No Mucosal Lesions)': -2, # Expects classic "volcano-like" ulcers
            'Inflammatory Labs (Leukocytosis/CRP)': 1
        }
    }

    print("Analyzing Clinical Case...\n")
    results = {}
    
    # Use pandas for a clean table display
    df_data = []

    for diagnosis, scores in scoring_matrix.items():
        total_score = 0
        calculation_str = []
        
        row = {'Diagnosis': diagnosis}
        for finding, finding_weight in findings.items():
            score_modifier = scores[finding]
            # The score for each finding is its presence (weight) times its relevance (modifier)
            # We simplify here by just using the modifier as the score contribution
            # since all findings are present.
            point = score_modifier
            total_score += point
            row[finding] = f"{point: d}" # Add finding score to the row
            # Format the string for the "equation"
            sign = '+' if point >= 0 else '-'
            calculation_str.append(f"({abs(point)})")
        
        results[diagnosis] = total_score
        row['Total Score'] = total_score
        df_data.append(row)
        
        print(f"Calculating score for: {diagnosis}")
        print(f"Equation: {' + '.join(calculation_str).replace('+ (-', '- (')}")
        print(f"Final Score = {total_score}\n")

    # Create and print the DataFrame
    df = pd.DataFrame(df_data)
    print("--- Scoring Summary ---")
    print(df.to_string(index=False))
    print("\n" + "="*25)

    # Determine the most likely diagnosis
    most_likely_diagnosis = max(results, key=results.get)
    print(f"Conclusion: The most likely diagnosis is '{most_likely_diagnosis}' with a score of {results[most_likely_diagnosis]}.")
    print("="*25)

if __name__ == '__main__':
    diagnose_esophageal_condition()