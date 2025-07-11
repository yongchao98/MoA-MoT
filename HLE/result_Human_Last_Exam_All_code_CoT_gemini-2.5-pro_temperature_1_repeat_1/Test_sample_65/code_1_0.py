def prioritize_patients():
    """
    Prioritizes spine surgery patients based on diagnosis and symptoms.

    The prioritization logic is as follows:
    1. Highest Priority (Score 3): Active neurologic deficits (e.g., cauda equina syndrome).
    2. High Priority (Score 2): Severe spinal instability with high risk of future neurologic injury (e.g., severe burst fracture).
    3. Lower Priority (Score 1): Less severe instability with no current neurologic deficits (e.g., compression fracture).
    """
    patients = [
        {
            "id": 1,
            "diagnosis": "Severe burst fracture of L2, no neurologic deficits",
            "priority_score": 2
        },
        {
            "id": 2,
            "diagnosis": "Compression fracture of L2 with mild traumatic spondylolisthesis of L1, no neurologic deficits",
            "priority_score": 1
        },
        {
            "id": 3,
            "diagnosis": "Split fracture of L2, with mildly disordered pelvic functions",
            "priority_score": 3
        }
    ]

    # Sort patients by priority_score in descending order
    sorted_patients = sorted(patients, key=lambda p: p['priority_score'], reverse=True)

    print("Surgical Priority from Highest to Lowest:")
    final_order = []
    for i, patient in enumerate(sorted_patients):
        print(f"{i+1}. Patient {patient['id']} (Score: {patient['priority_score']}) - {patient['diagnosis']}")
        final_order.append(f"Patient {patient['id']}")

    print("\nFinal Prioritization:")
    print(" -> ".join(final_order))

prioritize_patients()