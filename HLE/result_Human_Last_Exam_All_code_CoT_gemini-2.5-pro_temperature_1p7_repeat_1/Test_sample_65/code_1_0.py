def prioritize_patients():
    """
    Prioritizes spine surgery patients based on clinical urgency.
    Lower priority_rank means higher urgency (1 is the highest priority).
    """
    patients = [
        {
            "id": 1,
            "diagnosis": "Severe burst fracture of L2, no neurologic deficits",
            "reasoning": "High-grade spinal instability requires urgent stabilization.",
            "priority_rank": 2  # Second priority
        },
        {
            "id": 2,
            "diagnosis": "Compression fracture of L2 with mild traumatic spondylolisthesis of L1, no neurologic deficits",
            "reasoning": "Spinal instability, but less severe than a burst fracture.",
            "priority_rank": 3  # Third priority
        },
        {
            "id": 3,
            "diagnosis": "Split fracture of L2, with mildly disordered pelvic functions",
            "reasoning": "Active neurologic deficit (cauda equina syndrome) is a surgical emergency.",
            "priority_rank": 1  # Top priority
        }
    ]

    # Sort patients by their priority rank
    sorted_patients = sorted(patients, key=lambda p: p['priority_rank'])

    print("Patient Prioritization for Surgery (from highest to lowest):")
    final_order = []
    for patient in sorted_patients:
        print(f"\nPriority {patient['priority_rank']}: Patient {patient['id']}")
        print(f"  Diagnosis: {patient['diagnosis']}")
        print(f"  Reasoning: {patient['reasoning']}")
        final_order.append(str(patient['id']))
    
    print("\n-------------------------------------------------")
    print(f"Final Prioritization Order: Patient {final_order[0]}, Patient {final_order[1]}, Patient {final_order[2]}")
    print("-------------------------------------------------")


prioritize_patients()
<<<F>>>