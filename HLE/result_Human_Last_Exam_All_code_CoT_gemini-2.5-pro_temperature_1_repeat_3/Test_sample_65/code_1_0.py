def prioritize_spine_patients():
    """
    Prioritizes spine trauma patients based on clinical urgency.
    """
    # Patient data stored as a list of dictionaries
    patients = [
        {'id': 1, 'diagnosis': 'Severe burst fracture of L2', 'symptoms': 'no neurologic deficits'},
        {'id': 2, 'diagnosis': 'Compression fracture of L2 with mild traumatic spondylolisthesis of L1', 'symptoms': 'no neurologic deficits'},
        {'id': 3, 'diagnosis': 'Split fracture of L2', 'symptoms': 'mildly disordered pelvic functions'}
    ]

    def get_priority_score(patient):
        """
        Assigns a numerical score for sorting. Higher score means higher priority.
        - Priority 3: Active neurologic deficit (Cauda Equina Syndrome).
        - Priority 2: High mechanical instability with risk of neurologic injury.
        - Priority 1: Lower instability, no neurologic deficit.
        """
        symptoms = patient['symptoms']
        diagnosis = patient['diagnosis']

        if 'pelvic functions' in symptoms:
            return 3  # Neurosurgical emergency
        elif 'Severe burst fracture' in diagnosis:
            return 2  # Urgent due to high instability
        elif 'Compression fracture' in diagnosis:
            return 1  # Less urgent
        else:
            return 0 # Default

    # Sort the patients list based on the priority score in descending order
    sorted_patients = sorted(patients, key=get_priority_score, reverse=True)

    # Prepare the final output string showing the prioritization equation
    final_order = [p['id'] for p in sorted_patients]
    
    print("Surgical Priority from Highest to Lowest:")
    # The prompt requests to output each number in the final equation.
    # We will show the relationship as Patient 3 > Patient 1 > Patient 2
    print(f"Patient {final_order[0]} > Patient {final_order[1]} > Patient {final_order[2]}")

# Execute the function to print the result
prioritize_spine_patients()