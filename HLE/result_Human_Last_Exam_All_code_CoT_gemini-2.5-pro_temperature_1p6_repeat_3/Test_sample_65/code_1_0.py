def prioritize_patients():
    """
    This function prioritizes spine surgery patients based on their diagnosis
    and presence of neurologic deficits.
    """
    # Patient data with assigned priority scores.
    # Lower score means higher priority.
    # Priority 1: Active neurologic deficit (Cauda Equina Syndrome).
    # Priority 2: High-grade instability without neurologic deficit.
    # Priority 3: Lower-grade instability without neurologic deficit.
    patients = [
        {'id': 1, 'diagnosis': 'Severe burst fracture of L2, no neurologic deficits', 'priority': 2},
        {'id': 2, 'diagnosis': 'Compression fracture of L2 with mild traumatic spondylolisthesis of L1, no neurologic deficits', 'priority': 3},
        {'id': 3, 'diagnosis': 'Split fracture of L2, with mildly disordered pelvic functions', 'priority': 1}
    ]

    # Sort patients based on their priority score (ascending)
    sorted_patients = sorted(patients, key=lambda p: p['priority'])

    print("Surgical Priority List:")
    print("======================")

    final_order = []
    for i, patient in enumerate(sorted_patients):
        priority_label = f"{i+1}st"
        if i == 1: priority_label = "2nd"
        if i == 2: priority_label = "3rd"
        
        print(f"{priority_label} Priority:")
        print(f"  Patient {patient['id']}: {patient['diagnosis']}\n")
        final_order.append(str(patient['id']))
    
    # Per the instructions, output the numbers in the final equation/order
    print("Final Priority Order (Patient Numbers):")
    print(", ".join(final_order))

# Execute the function
prioritize_patients()