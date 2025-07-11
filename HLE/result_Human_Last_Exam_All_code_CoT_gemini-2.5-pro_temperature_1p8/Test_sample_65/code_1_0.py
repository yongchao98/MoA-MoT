def triage_spine_patients():
    """
    Prioritizes spine surgery patients based on clinical indications.
    The function assigns a priority score where a lower number indicates higher priority.
    - Priority 1: Active neurologic deficit (surgical emergency).
    - Priority 2: High degree of spinal instability with risk of future deficit.
    - Priority 3: Lower degree of spinal instability, no neurologic deficit.
    """

    patients = [
        {'id': 1, 'condition': 'Severe burst fracture of L2, no neurologic deficits', 'priority': 2},
        {'id': 2, 'condition': 'Compression fracture of L2 with mild traumatic spondylolisthesis, no neurologic deficits', 'priority': 3},
        {'id': 3, 'condition': 'Split fracture of L2, with mildly disordered pelvic functions', 'priority': 1}
    ]

    # Sort patients by their priority score (ascending)
    sorted_patients = sorted(patients, key=lambda x: x['priority'])

    print("Patient Triage Priority (from highest to lowest):")
    
    # Extract just the patient numbers for the final output
    prioritized_ids = [p['id'] for p in sorted_patients]
    
    # Print the final ordered list
    print(f"The final prioritization is: Patient {prioritized_ids[0]}, Patient {prioritized_ids[1]}, Patient {prioritized_ids[2]}")

# Run the triage
triage_spine_patients()