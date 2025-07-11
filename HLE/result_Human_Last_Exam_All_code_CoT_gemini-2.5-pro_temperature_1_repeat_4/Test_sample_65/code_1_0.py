def prioritize_patients():
    """
    Prioritizes three spine injury patients based on clinical urgency.
    """
    # Patient data with rationale for priority scoring
    # Lower score = higher priority
    patients = {
        'Patient 1': {
            'diagnosis': 'Severe burst fracture of L2, no neurologic deficits.',
            'priority_score': 2,
            'rationale': 'High degree of mechanical instability with risk of future neurologic injury.'
        },
        'Patient 2': {
            'diagnosis': 'Compression fracture of L2 with mild traumatic spondylolisthesis of L1, no neurologic deficits.',
            'priority_score': 3,
            'rationale': 'Least unstable fracture pattern with no neurologic deficits; lowest surgical urgency.'
        },
        'Patient 3': {
            'diagnosis': 'Split fracture of L2, with mildly disordered pelvic functions.',
            'priority_score': 1,
            'rationale': 'Active neurologic deficit (cauda equina syndrome), which is a surgical emergency.'
        }
    }

    # Sort patients based on their priority score (ascending)
    sorted_patients = sorted(patients.items(), key=lambda item: item[1]['priority_score'])

    print("Patient Prioritization for Surgery (Highest to Lowest):\n")

    for i, (patient_name, details) in enumerate(sorted_patients):
        # Extract the patient number from the name string
        patient_number = patient_name.split(' ')[1]
        print(f"Priority #{i+1}: Patient {patient_number}")
        print(f"  Diagnosis: {details['diagnosis']}")
        print(f"  Reasoning: {details['rationale']}\n")
    
    final_order = ", ".join([p_name for p_name, _ in sorted_patients])
    print(f"Final Triage Order: {final_order}")


prioritize_patients()