def triage_patients():
    """
    Prioritizes spine surgery patients based on clinical indications.
    """

    # Patient data with scoring for urgency.
    # Higher score means higher priority.
    # Scores are based on:
    # 1. Neurologic status (Cauda equina syndrome is highest priority)
    # 2. Fracture stability (More unstable fractures are higher priority)
    patients = [
        {
            "id": 1,
            "diagnosis": "Severe burst fracture of L2",
            "neuro_status": "No neurologic deficits",
            "urgency_score": 7  # High instability
        },
        {
            "id": 2,
            "diagnosis": "Compression fracture of L2 with mild traumatic spondylolisthesis of L1",
            "neuro_status": "No neurologic deficits",
            "urgency_score": 5  # Unstable, but less so than a severe burst fx
        },
        {
            "id": 3,
            "diagnosis": "Split fracture of L2",
            "neuro_status": "Mildly disordered pelvic functions",
            "urgency_score": 10 # Neurologic deficit (cauda equina) = surgical emergency
        }
    ]

    # Sort patients by urgency_score in descending order
    sorted_patients = sorted(patients, key=lambda x: x['urgency_score'], reverse=True)

    print("Prioritizing patients based on surgical indications (neurologic compromise and spinal instability):\n")

    priority_order = []
    for i, patient in enumerate(sorted_patients):
        print(f"Priority #{i+1}:")
        print(f"  Patient {patient['id']}")
        print(f"  Diagnosis: {patient['diagnosis']}")
        print(f"  Reasoning: {patient['neuro_status']}. Presence of neurologic deficit is the highest priority, followed by the degree of spinal instability.\n")
        priority_order.append(str(patient['id']))
        
    print("Final Prioritization Order (Top to Lowest):")
    # The final equation demonstrates the priority order.
    # Patient 3 is the top priority, followed by Patient 1, and then Patient 2.
    final_equation = " > ".join(priority_order)
    print(f"Patient {final_equation}")


triage_patients()
<<<F>>>