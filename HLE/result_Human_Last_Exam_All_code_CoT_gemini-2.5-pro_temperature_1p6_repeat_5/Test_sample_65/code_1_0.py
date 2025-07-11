def prioritize_spine_patients():
    """
    Prioritizes three spine injury patients based on clinical urgency.
    """
    # Patient data stored in a list of dictionaries
    patients = [
        {
            "id": 1,
            "description": "Severe burst fracture of L2, no neurologic deficits.",
            "reasoning": "High-grade mechanical instability creates a high risk of future neurologic injury. Urgent stabilization is needed."
        },
        {
            "id": 2,
            "description": "Compression fracture of L2 with mild traumatic spondylolisthesis of L1, no neurologic deficits.",
            "reasoning": "Indicates instability, but 'mild' and without deficits makes it less urgent than the others."
        },
        {
            "id": 3,
            "description": "Split fracture of L2, with mildly disordered pelvic functions.",
            "reasoning": "Active neurologic deficit (Cauda Equina Syndrome) is a surgical emergency requiring immediate intervention."
        }
    ]

    # Assign a priority score to each patient. A lower score indicates higher priority.
    for patient in patients:
        if "disordered pelvic functions" in patient["description"]:
            patient["priority_score"] = 1  # Top priority: Neurologic emergency
        elif "Severe burst fracture" in patient["description"]:
            patient["priority_score"] = 2  # High priority: High-grade instability
        else:
            patient["priority_score"] = 3  # Lower priority: Lower-grade instability

    # Sort patients based on their assigned priority score
    sorted_patients = sorted(patients, key=lambda x: x["priority_score"])

    # Print the final report
    print("Spine Patient Prioritization Report")
    print("====================================")
    for i, patient in enumerate(sorted_patients):
        print(f"Priority #{i+1}: Patient {patient['id']}")
        print(f"  - Condition: {patient['description']}")
        print(f"  - Rationale: {patient['reasoning']}")
        print("-" * 20)

    # Output the final prioritization as an "equation" with patient numbers
    final_order = " > ".join([str(p['id']) for p in sorted_patients])
    print(f"\nFinal Prioritization Order (Highest to Lowest): Patient {final_order}")

prioritize_spine_patients()