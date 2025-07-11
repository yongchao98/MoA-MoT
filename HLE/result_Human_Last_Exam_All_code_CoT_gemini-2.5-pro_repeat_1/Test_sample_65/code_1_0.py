def prioritize_patients():
    """
    This function prioritizes spine surgery patients based on clinical urgency.
    - Neurologic deficits are the highest priority.
    - Severe spinal instability is the second highest priority.
    """

    # Patient data with assigned priority scores
    # 3 = Highest priority (active neurologic deficit)
    # 2 = Medium priority (severe instability, high risk of deficit)
    # 1 = Lower priority (less severe instability, no deficit)
    patients = [
        {'id': 1, 'diagnosis': 'Severe burst fracture of L2, no neurologic deficits.', 'priority_score': 2},
        {'id': 2, 'diagnosis': 'Compression fracture of L2 with mild traumatic spondylolisthesis of L1, no neurologic deficits.', 'priority_score': 1},
        {'id': 3, 'diagnosis': 'Split fracture of L2, with mildly disordered pelvic functions.', 'priority_score': 3}
    ]

    # Sort patients by priority score in descending order
    sorted_patients = sorted(patients, key=lambda p: p['priority_score'], reverse=True)

    print("Prioritization of Spine Surgery Patients:")
    print("="*40)
    for i, patient in enumerate(sorted_patients, 1):
        print(f"Priority {i}: Patient {patient['id']}")
        print(f"  Diagnosis: {patient['diagnosis']}\n")
    
    # Extracting the numbers for the final equation-like output
    patient_ids = [p['id'] for p in sorted_patients]

    print("="*40)
    print("Final Priority Order Equation:")
    # The user requested to output each number in the final equation.
    # We will format this as Patient 3 > Patient 1 > Patient 2
    print(f"Patient {patient_ids[0]} > Patient {patient_ids[1]} > Patient {patient_ids[2]}")
    print("This corresponds to answer choice F.")

prioritize_patients()
<<<F>>>