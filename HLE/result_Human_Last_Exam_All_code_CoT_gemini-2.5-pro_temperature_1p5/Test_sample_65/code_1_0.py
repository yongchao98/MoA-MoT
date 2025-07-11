def triage_patients():
    """
    This function prioritizes spine surgery patients based on their diagnosis.

    Priority Levels:
    3 (Highest): Active neurologic deficit (e.g., cauda equina syndrome).
    2 (Medium): High-grade spinal instability without neurologic deficit (e.g., severe burst fracture).
    1 (Lowest): Moderate spinal instability without neurologic deficit (e.g., compression fracture with mild spondylolisthesis).
    """

    patients = [
        {'id': 1, 'diagnosis': 'Severe burst fracture of L2, no neurologic deficits.', 'priority_score': 2},
        {'id': 2, 'diagnosis': 'Compression fracture of L2 with mild traumatic spondylolisthesis of L1, no neurologic deficits.', 'priority_score': 1},
        {'id': 3, 'diagnosis': 'Split fracture of L2, with mildly disordered pelvic functions.', 'priority_score': 3}
    ]

    # Sort patients in descending order of their priority score
    sorted_patients = sorted(patients, key=lambda p: p['priority_score'], reverse=True)

    print("Surgical Prioritization from Highest to Lowest:")
    
    # Extract the patient IDs in the prioritized order
    prioritized_ids = [p['id'] for p in sorted_patients]
    
    # Print the final ordered list
    print(f"Patient {prioritized_ids[0]}, Patient {prioritized_ids[1]}, Patient {prioritized_ids[2]}")

triage_patients()
<<<F>>>