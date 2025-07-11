def triage_spine_patients():
    """
    Prioritizes three spine injury patients based on clinical urgency for surgery.
    """
    # Patient data with assigned priority scores. Lower score means higher priority.
    # Priority 1: Active neurological deficits (surgical emergency).
    # Priority 2: Highly unstable fracture with high risk of neurologic injury.
    # Priority 3: Less severe unstable fracture, no neurologic deficits.
    patients = [
        {
            "id": "Patient 1",
            "diagnosis": "Severe burst fracture of L2, no neurologic deficits.",
            "priority_score": 2
        },
        {
            "id": "Patient 2",
            "diagnosis": "Compression fracture of L2 with mild traumatic spondylolisthesis of L1, no neurologic deficits.",
            "priority_score": 3
        },
        {
            "id": "Patient 3",
            "diagnosis": "Split fracture of L2, with mildly disordered pelvic functions.",
            "priority_score": 1
        }
    ]

    # Sort patients based on their priority score (ascending)
    sorted_patients = sorted(patients, key=lambda x: x['priority_score'])

    print("Surgical Prioritization from Highest to Lowest Priority:")
    print("-------------------------------------------------------")
    for i, patient in enumerate(sorted_patients):
        print(f"Priority #{i+1}: {patient['id']}")
        print(f"   Diagnosis: {patient['diagnosis']}")
        print(f"   Reasoning: Patient {patient['id'].split(' ')[1]} is ranked #{i+1} due to the assigned priority score of {patient['priority_score']}.")
    print("-------------------------------------------------------")
    print(f"\nThe final prioritized order is: Patient {sorted_patients[0]['id'].split(' ')[1]}, Patient {sorted_patients[1]['id'].split(' ')[1]}, Patient {sorted_patients[2]['id'].split(' ')[1]}.")

triage_spine_patients()
<<<F>>>