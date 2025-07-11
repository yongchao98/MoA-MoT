def prioritize_spine_patients():
    """
    Analyzes and prioritizes three spine injury patients based on surgical indications.
    """
    # Define patients with attributes for prioritization.
    # neuro_priority: 1 = Emergency (pelvic dysfunction), 0 = No deficit.
    # stability_score: Higher score indicates greater instability.
    # 3: Severe Burst/Split, 2: Moderate, 1: Compression + Mild Spondylolisthesis.
    patients = [
        {
            "id": 1,
            "description": "Severe burst fracture of L2, no neurologic deficits.",
            "neuro_priority": 0,
            "stability_score": 3
        },
        {
            "id": 2,
            "description": "Compression fracture of L2 with mild traumatic spondylolisthesis of L1, no neurologic deficits.",
            "neuro_priority": 0,
            "stability_score": 1
        },
        {
            "id": 3,
            "description": "Split fracture of L2, with mildly disordered pelvic functions.",
            "neuro_priority": 1,
            "stability_score": 3 # A split fracture is also highly unstable.
        }
    ]

    # Sort patients. The primary key is neurological emergency status (descending).
    # The secondary key is the stability score (descending).
    # This places patients with neuro deficits first, then the most unstable ones.
    sorted_patients = sorted(
        patients,
        key=lambda p: (p['neuro_priority'], p['stability_score']),
        reverse=True
    )

    print("Surgical Prioritization Report")
    print("================================")
    print("Prioritization is based on: 1. Neurological Emergency, 2. Spinal Instability.\n")

    priority_order_string = ""
    for i, patient in enumerate(sorted_patients):
        print(f"Priority #{i+1}: Patient {patient['id']}")
        print(f"  - Diagnosis: {patient['description']}")
        if patient['neuro_priority'] == 1:
            print("  - Reasoning: Top priority due to active neurological compromise (pelvic dysfunction), indicating a surgical emergency.")
        else:
            if patient['stability_score'] == 3:
                print("  - Reasoning: High priority due to severe mechanical instability (burst/split fracture), which poses a risk of future neurologic injury.")
            else:
                print("  - Reasoning: Lower priority as the fracture pattern is the most stable among the cases, with no neurologic deficits.")
        print("-" * 20)
        
        # Build the final equation-like output string
        if i > 0:
            priority_order_string += " > "
        priority_order_string += f"Patient {patient['id']}"

    print(f"\nFinal Triage Order: {priority_order_string}")
    # The numbers in the final equation are: 3, 1, 2
    print(f"The numbers in the final prioritization are: {sorted_patients[0]['id']}, {sorted_patients[1]['id']}, {sorted_patients[2]['id']}")


prioritize_spine_patients()
<<<F>>>