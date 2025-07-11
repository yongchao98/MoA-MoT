def prioritize_patients():
    """
    Analyzes and prioritizes three spine injury patients based on surgical indications.
    """
    # Patient data with scores for neurologic deficit and instability.
    # Neurologic Deficit Score: 10 for present, 0 for absent.
    # Instability Score: 10 for severe, 5 for moderate/mild, 1 for stable.
    patients = [
        {
            "id": 1,
            "condition": "Severe burst fracture of L2, no neurologic deficits.",
            "neuro_score": 0,
            "instability_score": 10  # Severe instability
        },
        {
            "id": 2,
            "condition": "Compression fracture of L2 with mild traumatic spondylolisthesis of L1, no neurologic deficits.",
            "neuro_score": 0,
            "instability_score": 5   # Mild instability
        },
        {
            "id": 3,
            "condition": "Split fracture of L2, with mildly disordered pelvic functions.",
            "neuro_score": 10, # Active neurologic deficit (cauda equina syndrome)
            "instability_score": 8   # Split fractures are unstable
        }
    ]

    # Sort patients. The primary key is neurologic deficit, secondary is instability.
    # We sort in descending order to get the highest priority first.
    sorted_patients = sorted(patients, key=lambda p: (p['neuro_score'], p['instability_score']), reverse=True)

    print("Patient Prioritization for Surgery:")
    print("---------------------------------")
    print("Reasoning: Priority is given first to patients with active neurologic deficits, followed by those with severe spinal instability.")
    print("\nFinal Prioritization Order:")

    priority_list = []
    for i, patient in enumerate(sorted_patients):
        print(f"Priority {i+1}: Patient {patient['id']} ({patient['condition']})")
        priority_list.append(f"Patient {patient['id']}")
    
    # Print the final order in an equation-like format
    print("\nFinal Ranking:")
    print(" > ".join(priority_list))

prioritize_patients()
<<<F>>>