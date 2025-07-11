def prioritize_patients():
    """
    Prioritizes spine surgery patients based on neurologic status and fracture type.
    """
    # Define patients with their clinical information.
    # We assign scores for neurologic deficit and fracture severity.
    # Neurologic deficit is the highest priority (e.g., score of 10).
    # Fracture severity is ranked: severe burst > split > compression with mild slip.
    patients = [
        {
            'id': 1,
            'description': 'Severe burst fracture of L2, no neurologic deficits.',
            'neuro_score': 0,
            'instability_score': 8
        },
        {
            'id': 2,
            'description': 'Compression fracture of L2 with mild traumatic spondylolisthesis of L1, no neurologic deficits.',
            'neuro_score': 0,
            'instability_score': 4
        },
        {
            'id': 3,
            'description': 'Split fracture of L2, with mildly disordered pelvic functions.',
            'neuro_score': 10,  # Highest priority due to neurologic deficit
            'instability_score': 6
        }
    ]

    # Calculate total priority score for each patient
    for p in patients:
        p['total_score'] = p['neuro_score'] + p['instability_score']

    # Sort patients by total score in descending order (higher score = higher priority)
    sorted_patients = sorted(patients, key=lambda x: x['total_score'], reverse=True)

    print("Patient Prioritization for Surgery:\n")
    priority_order = []
    for i, p in enumerate(sorted_patients):
        print(f"Priority {i+1}:")
        print(f"  Patient {p['id']} - {p['description']}")
        print(f"  Reasoning: Neurologic Score = {p['neuro_score']}, Instability Score = {p['instability_score']}. Total Priority Score = {p['total_score']}.\n")
        priority_order.append(f"Patient {p['id']}")

    final_order = " -> ".join(priority_order)
    print("-----------------------------------------------------")
    print(f"Final Priority Order (Highest to Lowest): {final_order}")
    print("This corresponds to the sequence: Patient 3, Patient 1, Patient 2")


prioritize_patients()
<<<F>>>