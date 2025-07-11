import pandas as pd

def prioritize_patients():
    """
    This function prioritizes three spine injury patients for surgery based on clinical urgency.
    """
    # Patient data
    data = {
        'Patient ID': [1, 2, 3],
        'Diagnosis': [
            'Severe burst fracture of L2',
            'Compression fracture of L2 with mild traumatic spondylolisthesis of L1',
            'Split fracture of L2'
        ],
        'Neurologic Status': [
            'No neurologic deficits',
            'No neurologic deficits',
            'Mildly disordered pelvic functions'
        ]
    }
    patients_df = pd.DataFrame(data)

    # Triage logic: Higher score = higher priority
    # 3: Active neurologic deficit (Cauda Equina Syndrome risk)
    # 2: High mechanical instability, risk of future deficit
    # 1: Lower mechanical instability, no deficit
    priority_scores = []
    reasons = []
    for index, row in patients_df.iterrows():
        if "disordered pelvic functions" in row['Neurologic Status']:
            priority_scores.append(3)
            reasons.append("Highest priority due to active neurologic compromise (risk of Cauda Equina Syndrome).")
        elif "Severe burst fracture" in row['Diagnosis']:
            priority_scores.append(2)
            reasons.append("High priority due to severe spinal instability and risk of future neurologic deficit.")
        else:
            priority_scores.append(1)
            reasons.append("Lowest priority as the fracture is more stable with no neurologic deficits.")

    patients_df['Priority Score'] = priority_scores
    patients_df['Reasoning'] = reasons

    # Sort patients by priority
    sorted_patients = patients_df.sort_values(by='Priority Score', ascending=False)

    print("Patient Triage Priority Order:\n")
    # Using .to_string() to format the output nicely
    print(sorted_patients[['Patient ID', 'Priority Score', 'Reasoning']].to_string(index=False))
    
    # Extract the final order of patient numbers
    final_order = sorted_patients['Patient ID'].tolist()
    
    print("\n-------------------------------------------------")
    print("The final priority order is based on the patient IDs.")
    # The user requested to output each number in the final 'equation' or sequence
    print(f"Final Prioritized Sequence: Patient {final_order[0]}, Patient {final_order[1]}, Patient {final_order[2]}")
    print("-------------------------------------------------")


prioritize_patients()
<<<F>>>