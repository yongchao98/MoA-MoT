def prioritize_spine_patients():
    """
    Prioritizes spine surgery patients based on clinical urgency.
    The logic is based on standard neurosurgical and orthopedic triage principles:
    1. Active neurologic deficit (Cauda Equina Syndrome) is a surgical emergency.
    2. High-grade spinal instability (e.g., severe burst fracture) without deficit is urgent to prevent future injury.
    3. Lower-grade instability (e.g., compression fracture with mild listhesis) without deficit is the lowest priority among surgical candidates.
    """
    # Patient data with assigned priority scores (1=highest, 3=lowest)
    patients = [
        {
            "id": 1,
            "description": "Severe burst fracture of L2, no neurologic deficits.",
            "priority": 2,
            "reasoning": "High risk of future neurologic injury due to severe instability."
        },
        {
            "id": 2,
            "description": "Compression fracture of L2 with mild traumatic spondylolisthesis of L1, no neurologic deficits.",
            "priority": 3,
            "reasoning": "Instability present, but less severe than a burst fracture."
        },
        {
            "id": 3,
            "description": "Split fracture of L2, with mildly disordered pelvic functions.",
            "priority": 1,
            "reasoning": "Active neurologic deficit (Cauda Equina Syndrome), a surgical emergency."
        }
    ]

    # Sort patients based on their assigned priority
    sorted_patients = sorted(patients, key=lambda p: p['priority'])

    print("Surgical Prioritization (Highest to Lowest):")
    print("-" * 45)
    for patient in sorted_patients:
        print(f"Priority Rank {patient['priority']}: Patient {patient['id']}")
        print(f"  - Condition: {patient['description']}")
        print(f"  - Rationale: {patient['reasoning']}")
    
    # The final answer requires printing the numbers in the final order
    final_order = [p['id'] for p in sorted_patients]
    print("\nFinal Prioritization Order:")
    print(f"Patient {final_order[0]}, Patient {final_order[1]}, Patient {final_order[2]}")

prioritize_spine_patients()
<<<F>>>