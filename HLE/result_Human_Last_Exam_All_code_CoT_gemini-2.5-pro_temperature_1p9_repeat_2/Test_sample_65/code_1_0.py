def prioritize_patients():
    """
    Prioritizes three spine injury patients based on clinical urgency.
    """
    # Patient data with descriptions for clarity
    patients = {
        'Patient 1': {'diagnosis': 'Severe burst fracture of L2', 'neuro_status': 'no neurologic deficits'},
        'Patient 2': {'diagnosis': 'Compression fracture of L2 with mild traumatic spondylolisthesis of L1', 'neuro_status': 'no neurologic deficits'},
        'Patient 3': {'diagnosis': 'Split fracture of L2', 'neuro_status': 'mildly disordered pelvic functions'}
    }

    # Scoring system based on clinical principles
    # - Neurologic deficit is the highest priority (surgical emergency).
    # - Fracture instability determines priority in neurologically intact patients.
    neuro_scores = {
        'mildly disordered pelvic functions': 100,
        'no neurologic deficits': 0
    }

    instability_scores = {
        'Split fracture of L2': 50,  # Highly unstable 3-column injury
        'Severe burst fracture of L2': 40, # Highly unstable with risk to spinal canal
        'Compression fracture of L2 with mild traumatic spondylolisthesis of L1': 20 # Unstable, but less so than others
    }

    print("Calculating patient priority scores:\n")
    
    patient_priorities = []
    for name, data in patients.items():
        neuro_score = neuro_scores[data['neuro_status']]
        instability_score = instability_scores[data['diagnosis']]
        total_score = neuro_score + instability_score
        
        # Output the calculation for each patient
        print(f"{name} Priority Score = {neuro_score} (Neurologic Status) + {instability_score} (Fracture Instability) = {total_score}")
        
        patient_priorities.append({'name': name, 'score': total_score})

    # Sort patients by score in descending order (higher score = higher priority)
    sorted_patients = sorted(patient_priorities, key=lambda x: x['score'], reverse=True)

    print("\n--- Final Prioritization (from Highest to Lowest) ---")
    priority_order = []
    for i, patient in enumerate(sorted_patients):
        priority_order.append(patient['name'])
        print(f"{i+1}. {patient['name']} (Score: {patient['score']})")
    
    # This matches the reasoning: Patient 3, Patient 1, Patient 2
    final_order_string = ", ".join(priority_order)
    print(f"\nFinal Answer Order: {final_order_string}")


prioritize_patients()
<<<F>>>