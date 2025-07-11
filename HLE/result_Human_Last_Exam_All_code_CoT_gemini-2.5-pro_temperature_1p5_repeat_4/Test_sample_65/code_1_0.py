def prioritize_patients():
    """
    Prioritizes spine surgery patients based on diagnosis.
    
    The prioritization is based on standard clinical urgency:
    1. Highest Priority (Score 3): Active neurologic deficits (e.g., cauda equina syndrome).
    2. High Priority (Score 2): Severe spinal instability posing a high risk of neurologic injury.
    3. Lower Priority (Score 1): Less severe instability without current neurologic deficits.
    """
    
    patients = [
        {'id': 1, 'diagnosis': 'Severe burst fracture of L2, no neurologic deficits.'},
        {'id': 2, 'diagnosis': 'Compression fracture of L2 with mild traumatic spondylolisthesis of L1, no neurologic deficits.'},
        {'id': 3, 'diagnosis': 'Split fracture of L2, with mildly disordered pelvic functions.'}
    ]
    
    # Assign priority scores based on clinical keywords
    for patient in patients:
        if "pelvic functions" in patient['diagnosis']:
            patient['priority_score'] = 3
        elif "Severe burst fracture" in patient['diagnosis']:
            patient['priority_score'] = 2
        elif "Compression fracture" in patient['diagnosis']:
            patient['priority_score'] = 1
        else:
            patient['priority_score'] = 0 # Default/unknown

    # Sort patients by priority score in descending order
    sorted_patients = sorted(patients, key=lambda x: x['priority_score'], reverse=True)

    print("Spine Surgery Prioritization:")
    print("----------------------------")
    for i, patient in enumerate(sorted_patients):
        print(f"Priority #{i+1}: Patient {patient['id']} (Score: {patient['priority_score']}) - {patient['diagnosis']}")

    print("\nFinal Prioritization Order (Top to Lowest):")
    final_order = [p['id'] for p in sorted_patients]
    # The final output shows each number in the prioritized sequence
    print(f"Patient {final_order[0]}, Patient {final_order[1]}, Patient {final_order[2]}")

prioritize_patients()
<<<F>>>