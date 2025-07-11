def prioritize_patients():
    """
    This function prioritizes spine surgery patients based on clinical urgency.
    """
    # Patient data with scoring criteria:
    # Neurologic deficit (cauda equina) = 10 (Emergency)
    # Severe instability (burst fracture) = 7 (Urgent)
    # Moderate instability (compression + spondylolisthesis) = 4 (Less Urgent)
    patients = [
        {
            "id": 1,
            "diagnosis": "Severe burst fracture of L2",
            "neurology": "No neurologic deficits",
            "score": 7,
            "reason": "High mechanical instability without current nerve damage."
        },
        {
            "id": 2,
            "diagnosis": "Compression fracture of L2 with mild traumatic spondylolisthesis of L1",
            "neurology": "No neurologic deficits",
            "score": 4,
            "reason": "Less severe instability pattern and no nerve damage."
        },
        {
            "id": 3,
            "diagnosis": "Split fracture of L2",
            "neurology": "Mildly disordered pelvic functions",
            "score": 10,
            "reason": "Active neurologic compromise (cauda equina syndrome), a surgical emergency."
        }
    ]

    # Sort patients in descending order of their priority score
    prioritized_list = sorted(patients, key=lambda p: p['score'], reverse=True)

    print("Patient Triage Priority Order:\n")
    final_order = []
    for i, patient in enumerate(prioritized_list):
        print(f"Priority #{i+1}: Patient {patient['id']}")
        print(f"  - Score: {patient['score']}")
        print(f"  - Reason: {patient['reason']}")
        final_order.append(str(patient['id']))
    
    # Print the final order as a single line equation
    print("\nFinal Prioritization:")
    print(f"Patient {final_order[0]} > Patient {final_order[1]} > Patient {final_order[2]}")

prioritize_patients()