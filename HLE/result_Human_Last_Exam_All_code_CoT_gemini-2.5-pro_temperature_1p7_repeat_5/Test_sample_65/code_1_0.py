def prioritize_patients():
    """
    Prioritizes spine surgery patients based on clinical urgency.
    
    Scoring Logic:
    - Neurologic Deficit (pelvic dysfunction): +10 (Emergency)
    - Fracture Instability (Severe Burst): +5 (High risk)
    - Fracture Instability (Split): +4 (Unstable)
    - Fracture Instability (Compression/Mild Spondylo): +2 (Lower risk)
    """
    patients = [
        {'id': 1, 'diagnosis': 'Severe burst fracture of L2', 'neuro_status': 'no neurologic deficits'},
        {'id': 2, 'diagnosis': 'Compression fracture of L2 with mild traumatic spondylolisthesis of L1', 'neuro_status': 'no neurologic deficits'},
        {'id': 3, 'diagnosis': 'Split fracture of L2', 'neuro_status': 'mildly disordered pelvic functions'}
    ]

    def calculate_score(patient):
        score = 0
        # Score based on neurologic status
        if "pelvic functions" in patient['neuro_status']:
            score += 10
        
        # Score based on fracture type (instability)
        if "severe burst fracture" in patient['diagnosis']:
            score += 5
        elif "split fracture" in patient['diagnosis']:
            score += 4
        elif "compression fracture" in patient['diagnosis']:
            score += 2
            
        patient['priority_score'] = score
        return score

    # Sort patients by score in descending order
    patients.sort(key=calculate_score, reverse=True)

    print("Patient Prioritization from Highest to Lowest Urgency:")
    print("-----------------------------------------------------")
    final_order = []
    for patient in patients:
        print(f"Priority Rank #{len(final_order) + 1}: Patient {patient['id']} (Score: {patient['priority_score']})")
        print(f"  - Diagnosis: {patient['diagnosis']}")
        print(f"  - Neurology: {patient['neuro_status']}\n")
        final_order.append(f"Patient {patient['id']}")
    
    print("Final Prioritization Order:")
    print(f"{' > '.join(final_order)}")

prioritize_patients()
<<<F>>>