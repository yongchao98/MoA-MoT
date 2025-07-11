def prioritize_spine_patients():
    """
    This function prioritizes three spinal injury patients based on clinical urgency.
    """
    # Patient data with assigned priority scores.
    # Higher score = higher priority.
    # Scoring basis:
    # 10: Active neurologic deficit (surgical emergency)
    # 7: High mechanical instability without deficit (urgent)
    # 4: Moderate/mild instability without deficit (less urgent)
    patients = [
        {'id': 1, 'condition': 'Severe burst fracture, no neurologic deficits', 'score': 7},
        {'id': 2, 'condition': 'Compression fracture with mild spondylolisthesis, no neurologic deficits', 'score': 4},
        {'id': 3, 'condition': 'Split fracture, with mildly disordered pelvic functions', 'score': 10}
    ]

    # Sort patients in descending order based on their priority score
    sorted_patients = sorted(patients, key=lambda p: p['score'], reverse=True)

    print("Prioritizing patients based on surgical urgency (neurologic deficit > instability):\n")
    
    explanation = {
        3: "Highest Priority: Patient 3 has an active neurologic deficit (disordered pelvic functions), which is a surgical emergency (Cauda Equina Syndrome) requiring immediate intervention.",
        1: "Second Priority: Patient 1 has a severe burst fracture, indicating high mechanical instability. This is urgent to prevent future neurologic damage.",
        2: "Lowest Priority: Patient 2 has the most stable injury (compression fracture with mild slip) and no neurologic deficit, making this case the least urgent."
    }

    # Print the reasoning and the final ordered list
    final_order = []
    for patient in sorted_patients:
        patient_id = patient['id']
        print(explanation[patient_id])
        final_order.append(str(patient_id))
    
    print("\n---")
    print("The final prioritization from top to lowest priority is:")
    # The user requested to output each number in the final equation
    print(f"Patient {final_order[0]}, Patient {final_order[1]}, Patient {final_order[2]}")
    print("---\n")


prioritize_spine_patients()
# The final answer corresponds to the order: Patient 3, Patient 1, Patient 2
print("<<<F>>>")