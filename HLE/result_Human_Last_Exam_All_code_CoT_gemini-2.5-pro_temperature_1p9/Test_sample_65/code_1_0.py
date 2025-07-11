def prioritize_spine_patients():
    """
    Prioritizes spine surgery patients based on diagnosis and neurological status.
    The scoring system is as follows:
    - Neurologic Deficit: 100 points
    - No Neurologic Deficit: 0 points
    - Severe Burst Fracture: 50 points for instability
    - Split Fracture: 45 points for instability (a type of burst, also highly unstable)
    - Compression Fracture with Mild Spondylolisthesis: 20 points for instability
    """
    patients = [
        {'id': 1, 'diagnosis': 'Severe burst fracture', 'neurology': 'no neurologic deficits'},
        {'id': 2, 'diagnosis': 'Compression fracture with mild spondylolisthesis', 'neurology': 'no neurologic deficits'},
        {'id': 3, 'diagnosis': 'Split fracture', 'neurology': 'mildly disordered pelvic functions'}
    ]

    patient_scores = []

    print("Calculating patient priority scores...\n")

    for p in patients:
        neuro_score = 0
        instability_score = 0
        reasoning = []

        # Score based on neurology
        if "disordered" in p['neurology']:
            neuro_score = 100
            reasoning.append("Presence of neurologic deficit (100 points)")
        else:
            neuro_score = 0
            reasoning.append("No neurologic deficit (0 points)")

        # Score based on fracture instability
        if "Severe burst fracture" in p['diagnosis']:
            instability_score = 50
            reasoning.append("High instability from severe burst fracture (50 points)")
        elif "Split fracture" in p['diagnosis']:
            instability_score = 45
            reasoning.append("High instability from split fracture (45 points)")
        elif "Compression fracture" in p['diagnosis']:
            instability_score = 20
            reasoning.append("Moderate instability from compression fracture/mild spondylolisthesis (20 points)")
        
        total_score = neuro_score + instability_score
        
        print(f"Patient {p['id']}:")
        print(f"  Diagnosis: {p['diagnosis']}")
        print(f"  Neurology: {p['neurology']}")
        print(f"  Priority Score Calculation: {neuro_score} (neuro) + {instability_score} (instability) = {total_score}")
        print("-" * 20)
        
        p['score'] = total_score
        patient_scores.append(p)

    # Sort patients by score in descending order
    sorted_patients = sorted(patient_scores, key=lambda x: x['score'], reverse=True)

    print("\nFinal Prioritized List (Top to Lowest Priority):\n")
    for i, p in enumerate(sorted_patients, 1):
        print(f"{i}. Patient {p['id']} (Score: {p['score']})")

prioritize_spine_patients()
<<<F>>>