def prioritize_spine_patients():
    """
    Prioritizes three spine surgery patients based on clinical diagnosis.

    The prioritization is based on the urgency of surgical intervention:
    1.  Presence of neurologic deficits (especially cauda equina syndrome) is the highest priority.
    2.  High degree of spinal instability (e.g., severe burst fracture, split fracture) is the next priority.
    3.  Lesser degrees of instability without neurologic deficit are the lowest priority.
    """
    # Patient data with a scoring system for prioritization.
    # Higher score means higher priority.
    # Scores are based on: Neurologic status (100 for deficit) + Fracture Instability (scale 1-50)
    patients = [
        {
            "id": 1,
            "diagnosis": "Severe burst fracture of L2, no neurologic deficits.",
            "neuro_score": 0,
            "instability_score": 40  # Severe burst is highly unstable
        },
        {
            "id": 2,
            "diagnosis": "Compression fracture of L2 with mild traumatic spondylolisthesis of L1, no neurologic deficits.",
            "neuro_score": 0,
            "instability_score": 20  # Unstable, but described as 'mild'
        },
        {
            "id": 3,
            "diagnosis": "Split fracture of L2, with mildly disordered pelvic functions.",
            "neuro_score": 100, # Disordered pelvic functions = cauda equina signs = surgical emergency
            "instability_score": 50  # Split fracture is highly unstable
        }
    ]

    # Calculate total priority score for each patient
    for p in patients:
        p['priority_score'] = p['neuro_score'] + p['instability_score']

    # Sort patients by priority score in descending order
    sorted_patients = sorted(patients, key=lambda x: x['priority_score'], reverse=True)

    print("Patient Prioritization from Highest to Lowest Priority:")
    
    # Generate the final output string
    priority_order = []
    for patient in sorted_patients:
        priority_order.append(f"Patient {patient['id']}")
    
    print(" -> ".join(priority_order))
    
    print("\nReasoning:")
    print(f"1. Patient {sorted_patients[0]['id']}: Highest priority due to signs of cauda equina syndrome (disordered pelvic functions), which is a surgical emergency.")
    print(f"2. Patient {sorted_patients[1]['id']}: Second priority due to a severe burst fracture, indicating high mechanical instability and risk of future neurologic injury.")
    print(f"3. Patient {sorted_patients[2]['id']}: Lowest priority of the three, as the fracture and instability are described as 'mild' with no neurologic deficits.")

prioritize_spine_patients()
<<<F>>>