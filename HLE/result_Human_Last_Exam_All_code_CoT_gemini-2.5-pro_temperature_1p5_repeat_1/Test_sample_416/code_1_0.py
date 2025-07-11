import pandas as pd

def solve_diagnosis():
    """
    This function analyzes a clinical vignette by scoring potential diagnoses
    based on the evidence provided. It mimics a differential diagnosis process.
    """
    # Define the potential diagnoses
    diagnoses = {
        'A. Osteoarthritis': 0,
        'B. Charcot Arthropathy': 0,
        'C. Septic Arthritis': 0,
        'D. Chronic osteomyelitis': 0,
        'E. Pseudogout': 0
    }
    
    # Define the clinical findings and their corresponding scores for each diagnosis.
    # A positive score means the finding supports the diagnosis.
    # A negative score means the finding argues against the diagnosis.
    # A large negative score (-99) effectively rules out the diagnosis.
    clinical_findings = {
        'Onset after minor trauma (long walk)': {
            'A. Osteoarthritis': 1, 'B. Charcot Arthropathy': 2, 'C. Septic Arthritis': 0, 
            'D. Chronic osteomyelitis': 0, 'E. Pseudogout': 0, 'rationale': "Charcot is classically triggered by minor trauma."},
        'Acute inflammation (redness, swelling)': {
            'A. Osteoarthritis': 0, 'B. Charcot Arthropathy': 1, 'C. Septic Arthritis': 2, 
            'D. Chronic osteomyelitis': 1, 'E. Pseudogout': 2, 'rationale': "Consistent with inflammatory or infectious causes, and also Charcot."},
        'Negative X-rays initially': {
            'A. Osteoarthritis': -1, 'B. Charcot Arthropathy': 2, 'C. Septic Arthritis': 0, 
            'D. Chronic osteomyelitis': 0, 'E. Pseudogout': 0, 'rationale': "Very typical for early (Stage 0) Charcot; less so for established OA."},
        'Failure of NSAID treatment': {
            'A. Osteoarthritis': -1, 'B. Charcot Arthropathy': 2, 'C. Septic Arthritis': -1, 
            'D. Chronic osteomyelitis': -1, 'E. Pseudogout': -1, 'rationale': "Failure of anti-inflammatories points away from typical arthritis towards Charcot."},
        'Worsening on steroids': {
            'A. Osteoarthritis': -2, 'B. Charcot Arthropathy': 3, 'C. Septic Arthritis': -2, 
            'D. Chronic osteomyelitis': -2, 'E. Pseudogout': -2, 'rationale': "This is a key finding highly suggestive of Charcot arthropathy."},
        'Synovial fluid: No crystals': {
            'A. Osteoarthritis': 0, 'B. Charcot Arthropathy': 1, 'C. Septic Arthritis': 0, 
            'D. Chronic osteomyelitis': 0, 'E. Pseudogout': -99, 'rationale': "Effectively rules out Pseudogout (and Gout)."},
        'Synovial fluid: No organisms or WBCs': {
            'A. Osteoarthritis': 0, 'B. Charcot Arthropathy': 2, 'C. Septic Arthritis': -99, 
            'D. Chronic osteomyelitis': -1, 'E. Pseudogout': 0, 'rationale': "Rules out Septic Arthritis and strongly supports non-inflammatory Charcot."}
    }
    
    print("Evaluating clinical evidence to determine the most likely diagnosis...\n")
    
    score_calculations = {diag: [] for diag in diagnoses}

    for finding_name, scores in clinical_findings.items():
        print(f"Finding: {finding_name}")
        print(f"Rationale: {scores['rationale']}")
        for diag_name in diagnoses.keys():
            score = scores[diag_name]
            diagnoses[diag_name] += score
            score_calculations[diag_name].append(str(score))
        print("-" * 30)

    print("\nFinal Score Calculations:\n")
    for diag_name, score_list in score_calculations.items():
        equation = f" { ' + '.join(score_list) }".replace('+ -', '- ')
        print(f"{diag_name.ljust(25)} Score = {equation} = {diagnoses[diag_name]}")

    # Find the diagnosis with the highest score
    most_likely_diagnosis = max(diagnoses, key=diagnoses.get)

    print("\n-----------------------------------------------------------")
    print(f"Conclusion: The highest score belongs to '{most_likely_diagnosis}'.")
    print("This indicates Charcot Arthropathy is the most likely diagnosis based on the clinical evidence, particularly the failure of anti-inflammatory treatments and the bland synovial fluid analysis despite significant inflammation.")
    print("-----------------------------------------------------------")

solve_diagnosis()