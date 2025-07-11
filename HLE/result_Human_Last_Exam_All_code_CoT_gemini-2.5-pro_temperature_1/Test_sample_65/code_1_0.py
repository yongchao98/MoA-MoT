def prioritize_patients():
    """
    Prioritizes three spine injury patients based on surgical indications.

    Priority is determined by:
    1. Presence of neurologic deficits (highest priority).
    2. Degree of spinal instability (higher instability = higher priority).
    """

    # Patient data with assigned priority scores (1=highest, 3=lowest)
    # Patient 3: Highest priority due to active neurologic deficit (pelvic dysfunction).
    # Patient 1: Second priority due to severe mechanical instability (severe burst fracture).
    # Patient 2: Third priority due to instability that is less severe than Patient 1's.
    patients = [
        {'id': 1, 'condition': 'Severe burst fracture of L2, no neurologic deficits', 'priority': 2},
        {'id': 2, 'condition': 'Compression fracture of L2 with mild traumatic spondylolisthesis of L1, no neurologic deficits', 'priority': 3},
        {'id': 3, 'condition': 'Split fracture of L2, with mildly disordered pelvic functions', 'priority': 1}
    ]

    # Sort patients based on their priority score (ascending)
    sorted_patients = sorted(patients, key=lambda p: p['priority'])

    print("Patient Prioritization for Surgery (from highest to lowest priority):\n")

    final_order = []
    for patient in sorted_patients:
        final_order.append(f"Patient {patient['id']}")
        print(f"Priority {patient['priority']}: Patient {patient['id']}")
        print(f"  Diagnosis: {patient['condition']}")
        if patient['id'] == 3:
            print("  Rationale: Active neurologic deficit (pelvic dysfunction) indicates cauda equina involvement, a surgical emergency.\n")
        elif patient['id'] == 1:
            print("  Rationale: Severe mechanical instability from burst fracture poses a high risk of future neurologic compromise.\n")
        else:
            print("  Rationale: Unstable fracture, but less severe than a burst fracture and no active neurologic deficit.\n")
    
    # The final prioritized order is Patient 3, then Patient 1, then Patient 2.
    print("Final Prioritized Order:")
    print("Patient " + ", Patient ".join(map(str, [p['id'] for p in sorted_patients])))

prioritize_patients()