def solve_medical_case():
    """
    Analyzes a clinical vignette to determine the correct diagnosis.
    """
    # Step 1: Define the patient's data from the case description.
    patient_data = {
        'age': 27,
        'fever_duration_days': 4,
        'symptoms': [
            'fever', 'headaches', 'myalgia', 
            'disorientation', 'heart murmur'
        ],
        'history': 'camping trip to Oklahoma',
        'labs': {
            'Lyme IgM': 'elevated',
            'Lyme IgG': 'negative'
        }
    }

    # The organism that causes Lyme disease.
    lyme_causative_agent = "Borrelia burgdorferi"

    print("Analyzing the patient case:")
    print("---------------------------------")
    print(f"Patient Profile: A {patient_data['age']}-year-old with a {patient_data['fever_duration_days']}-day history of fever and other symptoms.")
    print(f"Key Findings: Neurological (disorientation) and cardiac (heart murmur) signs are present.")
    print(f"Crucial Lab Result: The Lyme serology shows an '{patient_data['labs']['Lyme IgM']}' IgM titer and a '{patient_data['labs']['Lyme IgG']}' IgG titer.")
    print("\nAnalysis:")
    print("An elevated IgM with a negative IgG is the classic serological profile for an acute (early) infection.")
    print(f"The disease associated with this lab finding is Lyme disease, which is caused by {lyme_causative_agent}.")
    
    print("\nDeriving the final answer:")
    
    # Final "equation" incorporating the numbers as requested.
    # It shows that the patient data plus the specific lab finding points directly to the causative agent.
    print(f"Patient Age ({patient_data['age']}) + Fever Duration ({patient_data['fever_duration_days']} days) + Lab Finding (Lyme IgM = '{patient_data['labs']['Lyme IgM']}') => Confirms an active infection with {lyme_causative_agent}.")

    print("\nConclusion:")
    print("The lab results explicitly state that the Lyme IgM titer is positive. Therefore, the titer for Borrelia burgdorferi is positive.")
    print("The correct answer choice is C.")

solve_medical_case()