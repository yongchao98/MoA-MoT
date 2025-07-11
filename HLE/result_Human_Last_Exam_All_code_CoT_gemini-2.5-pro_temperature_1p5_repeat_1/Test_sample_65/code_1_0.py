def triage_spine_patients():
    """
    Analyzes and prioritizes spine injury patients for surgery.
    """
    # Patient data with diagnoses
    patients = [
        {
            "id": 1,
            "diagnosis": "Severe burst fracture of L2, no neurologic deficits.",
            "reason": "High degree of spinal instability from a severe burst fracture creates a significant risk of future neurologic injury."
        },
        {
            "id": 2,
            "diagnosis": "Compression fracture of L2 with mild traumatic spondylolisthesis of L1, no neurologic deficits.",
            "reason": "The fracture pattern and mild instability without neurologic deficits represent the least urgent scenario."
        },
        {
            "id": 3,
            "diagnosis": "Split fracture of L2, with mildly disordered pelvic functions.",
            "reason": "Active neurologic deficit (disordered pelvic functions) is a surgical emergency (cauda equina syndrome)."
        }
    ]

    # Assign priority scores: 1 (highest), 2 (medium), 3 (lowest)
    for patient in patients:
        if "pelvic functions" in patient["diagnosis"]:
            patient["priority_score"] = 1  # Emergency
        elif "Severe burst fracture" in patient["diagnosis"]:
            patient["priority_score"] = 2  # Urgent due to instability
        else:
            patient["priority_score"] = 3  # Less urgent

    # Sort patients by priority score (ascending)
    sorted_patients = sorted(patients, key=lambda p: p["priority_score"])

    # Print the results
    print("Spine Surgery Triage Prioritization:")
    print("=" * 40)
    for i, patient in enumerate(sorted_patients):
        print(f"Priority #{i+1}: Patient {patient['id']}")
        print(f"  Diagnosis: {patient['diagnosis']}")
        print(f"  Reasoning: {patient['reason']}")
        print("-" * 40)
    
    final_order = [p['id'] for p in sorted_patients]
    print(f"The final priority order is Patient {final_order[0]}, then Patient {final_order[1]}, then Patient {final_order[2]}.")

# Execute the function
triage_spine_patients()
