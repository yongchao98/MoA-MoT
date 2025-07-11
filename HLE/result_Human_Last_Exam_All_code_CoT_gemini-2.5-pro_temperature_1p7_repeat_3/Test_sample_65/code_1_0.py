def triage_spine_patients():
    """
    Prioritizes spine injury patients based on surgical indications.
    The function creates patient profiles, assigns a priority score,
    and then prints the reasoned prioritization order.
    """
    
    # Patient data stored in a list of dictionaries
    patients = [
        {
            'id': 1,
            'diagnosis': 'Severe burst fracture of L2, no neurologic deficits.',
            'priority_score': 2, # Second priority due to high mechanical instability
            'reasoning': 'A severe burst fracture indicates high mechanical instability. There is a significant risk of future neurologic injury if the spine is not stabilized. This is urgent.'
        },
        {
            'id': 2,
            'diagnosis': 'Compression fracture of L2 with mild traumatic spondylolisthesis of L1, no neurologic deficits.',
            'priority_score': 3, # Third priority due to instability, but less severe than a burst fracture
            'reasoning': 'The combination of fracture and spondylolisthesis indicates instability requiring surgery, but it is less acutely dangerous than a severe burst fracture or a case with active neurologic deficits.'
        },
        {
            'id': 3,
            'diagnosis': 'Split fracture of L2, with mildly disordered pelvic functions.',
            'priority_score': 1, # Top priority due to neurologic deficit
            'reasoning': "'Disordered pelvic functions' signifies cauda equina syndrome, a neurosurgical emergency requiring immediate decompression to prevent permanent nerve damage."
        }
    ]

    # Sort patients based on their priority score (1 is highest priority)
    sorted_patients = sorted(patients, key=lambda p: p['priority_score'])

    print("Spine Patient Surgical Prioritization:\n")
    
    # Print the detailed reasoning for each patient in order of priority
    for i, patient in enumerate(sorted_patients):
        print(f"Priority #{i+1}: Patient {patient['id']}")
        print(f"  - Diagnosis: {patient['diagnosis']}")
        print(f"  - Reasoning: {patient['reasoning']}\n")

    # Print the final equation showing the order of patient numbers
    final_order = [str(p['id']) for p in sorted_patients]
    print("Final Prioritization from top to lowest priority:")
    print("Patient " + ", Patient ".join(final_order))

# Run the triage function
triage_spine_patients()