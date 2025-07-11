def prioritize_patients():
    """
    Prioritizes spine surgery patients based on clinical urgency.

    The priority is determined by the following hierarchy:
    1. Presence of neurological deficits (Highest Priority).
    2. High-degree of spinal instability.
    3. Lower-degree of spinal instability.
    """
    # Patient data with their diagnosis and priority score (1=highest, 3=lowest)
    patients = [
        {
            "id": 1,
            "condition": "Severe burst fracture of L2, no neurologic deficits",
            "priority_score": 2, # High instability, but no neuro deficit
            "reasoning": "High-risk mechanical instability requires urgent stabilization to prevent future neurological injury."
        },
        {
            "id": 2,
            "condition": "Compression fracture of L2 with mild traumatic spondylolisthesis of L1, no neurologic deficits",
            "priority_score": 3, # Less severe instability
            "reasoning": "The least urgent case; instability is present but described as 'mild' with no neurological symptoms."
        },
        {
            "id": 3,
            "condition": "Split fracture of L2, with mildly disordered pelvic functions",
            "priority_score": 1, # Active neurological deficit
            "reasoning": "Active neurologic deficit (cauda equina syndrome) is a surgical emergency requiring immediate decompression."
        }
    ]

    # Sort patients based on their priority score (ascending)
    prioritized_list = sorted(patients, key=lambda p: p['priority_score'])

    print("Patient Prioritization for Surgery (Highest to Lowest):")
    print("-" * 55)
    for i, patient in enumerate(prioritized_list, 1):
        print(f"Priority #{i}: Patient {patient['id']}")
        print(f"  Condition: {patient['condition']}")
        print(f"  Reasoning: {patient['reasoning']}\n")
    
    final_order = ", ".join([f"Patient {p['id']}" for p in prioritized_list])
    print(f"Final Prioritized Order: {final_order}")


prioritize_patients()
