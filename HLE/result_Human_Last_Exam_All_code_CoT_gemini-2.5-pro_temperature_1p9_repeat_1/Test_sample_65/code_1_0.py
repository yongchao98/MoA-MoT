def prioritize_patients():
    """
    This function prioritizes spine surgery patients based on their clinical presentation.
    It assigns a priority score and then sorts the patients.
    """
    # Define patient data. Priority scores are assigned based on clinical urgency:
    # - Neurologic deficit (especially pelvic function) is the highest priority.
    # - Severe mechanical instability (burst fracture) is the second highest.
    # - Milder instability is the lowest priority.
    patients = [
        {
            'id': 1,
            'diagnosis': 'Severe burst fracture of L2, no neurologic deficits.',
            'priority_score': 50  # High score for severe instability
        },
        {
            'id': 2,
            'diagnosis': 'Compression fracture of L2 with mild traumatic spondylolisthesis of L1, no neurologic deficits.',
            'priority_score': 20  # Lower score for less severe instability
        },
        {
            'id': 3,
            'diagnosis': 'Split fracture of L2, with mildly disordered pelvic functions.',
            'priority_score': 100 # Highest score for active neurologic deficit
        }
    ]

    # Sort patients by priority_score in descending order
    sorted_patients = sorted(patients, key=lambda x: x['priority_score'], reverse=True)

    print("Patient Prioritization for Surgery:\n")
    print("Reasoning: The primary factor for surgical urgency is the presence of neurologic deficits. The secondary factor is the degree of spinal instability.\n")
    
    final_order = []
    for i, patient in enumerate(sorted_patients):
        priority_level = ""
        if i == 0:
            priority_level = "Top Priority"
        elif i == 1:
            priority_level = "Second Priority"
        else:
            priority_level = "Lowest Priority"
            
        print(f"{priority_level}: Patient {patient['id']}")
        print(f"  - Diagnosis: {patient['diagnosis']}")
        if patient['id'] == 3:
            print("  - Justification: Presence of neurologic deficit (disordered pelvic functions) is a surgical emergency.\n")
        elif patient['id'] == 1:
            print("  - Justification: Severe mechanical instability from a burst fracture poses a high risk of future neurologic injury and deformity.\n")
        else:
            print("  - Justification: Lower degree of instability and no neurologic deficits make this case less urgent.\n")
        
        final_order.append(f"Patient {patient['id']}")
        
    print("---------------------------------")
    print("Final Priority Order: " + ", ".join(final_order))
    print("This corresponds to the answer choice: F")

prioritize_patients()
<<<F>>>