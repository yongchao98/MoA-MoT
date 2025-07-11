def prioritize_patients():
    """
    Prioritizes spine surgery patients based on clinical urgency.

    Priority is determined by:
    1. Presence of Neurologic Deficits (Highest Priority)
    2. Degree of Spinal Instability (Secondary Factor)
    """

    # Patient data with descriptions and assigned priority scores
    # A lower score means higher priority.
    # Score 1: Active neurologic deficit (surgical emergency)
    # Score 2: High-grade instability without neurologic deficit
    # Score 3: Lower-grade instability without neurologic deficit
    patients = [
        {'id': 1, 'diagnosis': 'Severe burst fracture of L2, no neurologic deficits.', 'priority_score': 2},
        {'id': 2, 'diagnosis': 'Compression fracture of L2 with mild traumatic spondylolisthesis of L1, no neurologic deficits.', 'priority_score': 3},
        {'id': 3, 'diagnosis': 'Split fracture of L2, with mildly disordered pelvic functions.', 'priority_score': 1}
    ]

    # Sort patients based on their priority score (ascending)
    sorted_patients = sorted(patients, key=lambda p: p['priority_score'])

    print("Spine Surgery Triage - Prioritization List:")
    print("-" * 45)
    print("Priority | Patient | Diagnosis")
    print("-" * 45)

    for i, patient in enumerate(sorted_patients):
        priority_rank = i + 1
        print(f"{priority_rank:<8} | Patient {patient['id']} | {patient['diagnosis']}")
    
    final_order = ", ".join([f"Patient {p['id']}" for p in sorted_patients])
    print("\nFinal Prioritization from Highest to Lowest:")
    print(final_order)

prioritize_patients()