def prioritize_patients():
    """
    Prioritizes spine surgery patients based on diagnosis and neurological status.

    Priority is determined by:
    1. Presence of neurological deficits (highest priority).
    2. Degree of spinal instability (severe burst fracture > compression fx with listhesis).
    """

    # Patient data with assigned priority scores based on clinical urgency
    # Score 1: Neurological deficit (surgical emergency)
    # Score 2: High mechanical instability without deficits
    # Score 3: Moderate mechanical instability without deficits
    patients = [
        {'id': 1, 'condition': 'Severe burst fracture of L2, no neurologic deficits', 'priority_score': 2},
        {'id': 2, 'condition': 'Compression fracture of L2 with mild traumatic spondylolisthesis, no neurologic deficits', 'priority_score': 3},
        {'id': 3, 'condition': 'Split fracture of L2, with mildly disordered pelvic functions', 'priority_score': 1},
    ]

    # Sort patients by their priority score (ascending)
    sorted_patients = sorted(patients, key=lambda p: p['priority_score'])

    # Print the final prioritization list
    print("Surgical Prioritization from Highest to Lowest:")
    priority_order = []
    for patient in sorted_patients:
        priority_order.append(str(patient['id']))
    
    # The final equation is the ordered list of patients
    print("Patient " + ", Patient ".join(priority_order))

prioritize_patients()
<<<F>>>