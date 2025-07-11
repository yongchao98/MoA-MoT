def prioritize_patients():
    """
    Prioritizes spine surgery patients based on clinical diagnosis.

    Priority is determined by:
    1. Presence of neurologic deficits (highest priority).
    2. Degree of spinal instability (higher instability is higher priority).
    """

    # Patient data with assigned priority scores (1 = highest priority)
    # Patient 3: Neurologic deficit -> Priority 1
    # Patient 1: Severe instability (burst fracture) -> Priority 2
    # Patient 2: Moderate instability (compression fx + mild spondylolisthesis) -> Priority 3
    patients = [
        {'id': 1, 'diagnosis': 'Severe burst fracture of L2, no neurologic deficits', 'priority_score': 2},
        {'id': 2, 'diagnosis': 'Compression fracture of L2 with mild traumatic spondylolisthesis of L1, no neurologic deficits', 'priority_score': 3},
        {'id': 3, 'diagnosis': 'Split fracture of L2, with mildly disordered pelvic functions', 'priority_score': 1}
    ]

    # Sort patients based on their priority score (ascending)
    sorted_patients = sorted(patients, key=lambda p: p['priority_score'])

    print("Patient Prioritization for Surgery:")
    print("---------------------------------")
    
    # We will build the final equation string as requested
    final_order_list = []
    
    for i, patient in enumerate(sorted_patients):
        print(f"Priority #{i+1}: Patient {patient['id']}")
        print(f"  Diagnosis: {patient['diagnosis']}")
        if patient['id'] == 3:
            print("  Reason: Highest priority due to active neurologic deficits (cauda equina syndrome), which is a surgical emergency.")
        elif patient['id'] == 1:
            print("  Reason: High priority due to severe spinal instability from burst fracture, posing a risk of future neurologic injury and deformity.")
        else:
            print("  Reason: Lower priority. While unstable, the degree of instability is less severe than Patient 1.")
        
        final_order_list.append(f"Patient {patient['id']}")

    # Print the final prioritized order as an "equation"
    print("\nFinal Prioritization Order:")
    print(" > ".join(final_order_list))

prioritize_patients()