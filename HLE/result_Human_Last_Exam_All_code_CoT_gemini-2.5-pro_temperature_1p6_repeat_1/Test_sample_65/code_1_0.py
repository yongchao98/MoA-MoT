def triage_spine_patients():
    """
    Prioritizes three spine injury patients based on surgical indications.
    The function prints the reasoning and the final prioritized order.
    """

    # Patient data with key findings for prioritization
    patients = [
        {'id': 1, 'diagnosis': 'Severe burst fracture of L2', 'neuro_status': 'no neurologic deficits'},
        {'id': 2, 'diagnosis': 'Compression fracture of L2 with mild traumatic spondylolisthesis of L1', 'neuro_status': 'no neurologic deficits'},
        {'id': 3, 'diagnosis': 'Split fracture of L2', 'neuro_status': 'mildly disordered pelvic functions'}
    ]

    print("Spine Surgery Triage Protocol")
    print("--------------------------------")
    print("Prioritization Criteria:")
    print("1. Presence of Neurologic Deficit (especially cauda equina syndrome) -> Highest Priority")
    print("2. Degree of Mechanical Instability -> Second Highest Priority")
    print("--------------------------------\n")

    # Assign a priority score to each patient
    # Higher score means higher priority
    for patient in patients:
        if "disordered pelvic functions" in patient['neuro_status']:
            # Active cauda equina syndrome is a surgical emergency
            patient['priority_score'] = 3
            patient['reasoning'] = "Active neurologic deficit (cauda equina syndrome) requires urgent surgical decompression."
        elif "Severe burst fracture" in patient['diagnosis']:
            # High degree of mechanical instability poses risk of future neurologic injury
            patient['priority_score'] = 2
            patient['reasoning'] = "Severe mechanical instability from a burst fracture requires stabilization to prevent future decline."
        else:
            # Mild instability without neurologic deficits is the lowest priority among these cases
            patient['priority_score'] = 1
            patient['reasoning'] = "Relatively stable fracture pattern with no neurologic deficits; least urgent case."

    # Sort patients by priority score in descending order
    sorted_patients = sorted(patients, key=lambda p: p['priority_score'], reverse=True)

    print("Triage Result (from top priority to lowest priority):\n")
    final_order = []
    for i, patient in enumerate(sorted_patients):
        print(f"Priority #{i+1}: Patient {patient['id']}")
        print(f"  Diagnosis: {patient['diagnosis']}")
        print(f"  Clinical Finding: {patient['neuro_status']}")
        print(f"  Reasoning: {patient['reasoning']}\n")
        final_order.append(str(patient['id']))
        
    print("--------------------------------")
    print(f"Final Prioritization Order: Patient {final_order[0]}, Patient {final_order[1]}, Patient {final_order[2]}")
    
    # Matching the order (3, 1, 2) to the provided answer choices
    answer_choice = "F"
    print(f"This corresponds to Answer Choice: {answer_choice}")


# Execute the triage process
triage_spine_patients()

# Return the final answer in the specified format
print("<<<F>>>")
