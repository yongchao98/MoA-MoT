def prioritize_spine_patients():
    """
    Prioritizes three spine injury patients based on clinical urgency
    using a scoring system that reflects neurosurgical principles.
    """
    # 1. Define patient data in a structured format
    patients = [
        {
            'id': 1,
            'diagnosis': 'Severe burst fracture of L2',
            'neuro_status': 'no neurologic deficits'
        },
        {
            'id': 2,
            'diagnosis': 'Compression fracture of L2 with mild traumatic spondylolisthesis of L1',
            'neuro_status': 'no neurologic deficits'
        },
        {
            'id': 3,
            'diagnosis': 'Split fracture of L2',
            'neuro_status': 'mildly disordered pelvic functions'
        }
    ]

    # 2. Define the scoring logic in a function
    def calculate_priority_score(patient):
        """Calculates a priority score based on neurologic status and fracture type."""
        # Neurologic score: Cauda equina signs are a surgical emergency.
        neuro_score = 100 if "pelvic functions" in patient['neuro_status'] else 0

        # Instability score: Based on fracture pattern.
        # Split fractures (3-column) are highly unstable.
        # Severe burst fractures are highly unstable.
        # Compression fractures are generally less unstable.
        instability_score = 0
        if "Split fracture" in patient['diagnosis']:
            instability_score = 60
        elif "Severe burst fracture" in patient['diagnosis']:
            instability_score = 50
        elif "Compression fracture" in patient['diagnosis']:
            instability_score = 20
        
        # Total score is the sum of both components
        total_score = neuro_score + instability_score
        
        return {
            'neuro': neuro_score,
            'instability': instability_score,
            'total': total_score
        }

    # 3. Calculate scores for each patient and store them
    for patient in patients:
        scores = calculate_priority_score(patient)
        patient['scores'] = scores

    # 4. Sort patients based on total score in descending order
    sorted_patients = sorted(patients, key=lambda p: p['scores']['total'], reverse=True)

    # 5. Print the detailed results and final conclusion
    print("--- Patient Prioritization Analysis ---\n")
    print("Urgency is determined by scoring neurologic deficits and spinal instability.")
    print("Equation: Priority Score = Neurologic Score + Instability Score\n")

    print("--- Calculated Patient Scores ---")
    for p in patients:
        print(f"Patient {p['id']}: Neurologic Score ({p['scores']['neuro']}) + Instability Score ({p['scores']['instability']}) = Total Priority Score: {p['scores']['total']}")

    print("\n--- Final Prioritization Order ---")
    rank = 1
    final_order_list = []
    for p in sorted_patients:
        final_order_list.append(f"Patient {p['id']}")
        print(f"{rank}. Patient {p['id']} (Score: {p['scores']['total']}) - Reason: {p['diagnosis']} with {p['neuro_status']}.")
        rank += 1
    
    final_order_string = ", ".join(final_order_list)
    print(f"\nThe determined priority order is: {final_order_string}.")
    print("This corresponds to answer choice F.")

# Execute the function to get the answer
prioritize_spine_patients()