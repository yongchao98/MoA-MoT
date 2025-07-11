def prioritize_spine_patients():
    """
    Analyzes and prioritizes three spine injury patients based on surgical urgency.
    """
    # Patient data representation: [ID, Diagnosis, Neurologic Status, Priority Score]
    # Priority Score: 3 = Highest, 2 = Medium, 1 = Lowest
    patients = [
        {"id": 1, "diagnosis": "Severe burst fracture of L2", "neuro_status": "No neurologic deficits", "priority_score": 2, "reasoning": "Severe mechanical instability poses a high risk of future neurologic injury, requiring urgent stabilization."},
        {"id": 2, "diagnosis": "Compression fracture of L2 with mild traumatic spondylolisthesis of L1", "neuro_status": "No neurologic deficits", "priority_score": 1, "reasoning": "Lowest degree of instability and no neurologic deficit; least urgent need for surgery."},
        {"id": 3, "diagnosis": "Split fracture of L2", "neuro_status": "Mildly disordered pelvic functions", "priority_score": 3, "reasoning": "Presence of a neurologic deficit (cauda equina syndrome) is a surgical emergency requiring immediate attention."}
    ]

    # Sort patients by priority score in descending order
    sorted_patients = sorted(patients, key=lambda x: x['priority_score'], reverse=True)

    print("Spine Patient Triage Prioritization:")
    print("------------------------------------\n")
    print("Reasoning: Prioritization is based on the presence of active neurologic deficits, followed by the degree of spinal instability.\n")

    priority_number = 1
    final_order = []
    for patient in sorted_patients:
        print(f"Priority #{priority_number}: Patient {patient['id']}")
        print(f"  - Diagnosis: {patient['diagnosis']}")
        print(f"  - Neurologic Status: {patient['neuro_status']}")
        print(f"  - Justification: {patient['reasoning']}\n")
        final_order.append(f"Patient {patient['id']}")
        priority_number += 1
    
    print("------------------------------------")
    print(f"Final Prioritization Order (Highest to Lowest):")
    print(" -> ".join(final_order))
    
prioritize_spine_patients()
<<<F>>>