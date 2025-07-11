import pandas as pd

def diagnose_hypoxemia_cause():
    """
    Analyzes a clinical vignette to determine the most likely cause of hypoxemia
    by scoring diagnoses against key patient findings.
    """

    # Key clinical findings from the vignette
    patient_findings = {
        'presentation': 'ARDS-like (severe hypoxemia + bilateral crackles)',
        'timeline': 'subacute (29 days post-op)',
        'context': 'major abdominal surgery (Whipple)'
    }

    # Potential diagnoses with characteristics for scoring
    # Scoring criteria:
    # - fits_presentation: Does it cause ARDS? (+5)
    # - fits_timeline: Is the timeline plausible? (+5 for subacute, -5 for acute)
    # - fits_context: Is it a known complication of major surgery? (+5)
    diagnoses_data = {
        'A': {'Diagnosis': 'Acute blood transfusion reaction', 'fits_presentation': 5, 'fits_timeline': -5, 'fits_context': 0},
        'B': {'Diagnosis': 'Iodine-related reaction', 'fits_presentation': 0, 'fits_timeline': -5, 'fits_context': 0},
        'C': {'Diagnosis': 'Sensitivity reaction', 'fits_presentation': 1, 'fits_timeline': 0, 'fits_context': 0},
        'D': {'Diagnosis': 'Sepsis', 'fits_presentation': 5, 'fits_timeline': 5, 'fits_context': 5},
        'E': {'Diagnosis': 'Myocyte necrosis', 'fits_presentation': 1, 'fits_timeline': 0, 'fits_context': 1},
        'F': {'Diagnosis': 'Respiratory deconditioning', 'fits_presentation': 0, 'fits_timeline': 2, 'fits_context': 2},
        'G': {'Diagnosis': 'Lung exhaustion', 'fits_presentation': 0, 'fits_timeline': 0, 'fits_context': 0},
        'H': {'Diagnosis': 'Air pollution sensitivity', 'fits_presentation': 0, 'fits_timeline': 0, 'fits_context': -5}
    }

    results = []
    print("Evaluating diagnoses based on patient's clinical picture...\n")

    for key, data in diagnoses_data.items():
        score = data['fits_presentation'] + data['fits_timeline'] + data['fits_context']
        equation = f"{data['fits_presentation']} (presentation) + {data['fits_timeline']} (timeline) + {data['fits_context']} (context)"
        results.append({
            'Option': key,
            'Diagnosis': data['Diagnosis'],
            'Likelihood Score': score,
            'Scoring Equation': equation
        })

    # Create a DataFrame for nice printing
    df = pd.DataFrame(results).set_index('Option')
    
    # Print the equation and result for each diagnosis
    print("Likelihood Calculation (Score = Presentation Fit + Timeline Fit + Context Fit):")
    for index, row in df.iterrows():
      print(f"Option {index} ({row['Diagnosis']}):")
      print(f"  Final Score = {row['Scoring Equation']} = {row['Likelihood Score']}")
      print("-" * 30)

    # Determine the most likely diagnosis
    most_likely = df['Likelihood Score'].idxmax()
    
    print(f"\nConclusion: The diagnosis with the highest likelihood score is Option {most_likely}.")
    print(f"The clinical picture is most consistent with Sepsis leading to ARDS.")

    return most_likely

if __name__ == '__main__':
    final_answer = diagnose_hypoxemia_cause()
    # The final answer is wrapped in the required format.
    print(f"\n<<<{final_answer}>>>")
