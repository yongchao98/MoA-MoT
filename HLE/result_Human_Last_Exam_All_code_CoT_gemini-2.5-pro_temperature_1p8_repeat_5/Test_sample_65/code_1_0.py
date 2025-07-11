import sys

def solve_triage_task():
    """
    This function prioritizes three spine surgery patients based on their clinical presentation.
    """
    # Define patient data
    patients = [
        {
            "id": 1,
            "diagnosis": "Severe burst fracture of L2, no neurologic deficits.",
            "neuro_deficit": False,
            "instability": "Severe"
        },
        {
            "id": 2,
            "diagnosis": "Compression fracture of L2 with mild traumatic spondylolisthesis of L1, no neurologic deficits.",
            "neuro_deficit": False,
            "instability": "Moderate" # Spondylolisthesis implies instability, but it's less severe than a "severe burst"
        },
        {
            "id": 3,
            "diagnosis": "Split fracture of L2, with mildly disordered pelvic functions.",
            "neuro_deficit": True, # This is a sign of Cauda Equina Syndrome, a surgical emergency.
            "instability": "High"
        }
    ]

    # Assign a priority score based on clinical urgency
    for patient in patients:
        if patient["neuro_deficit"]:
            # Active neurologic deficit is the highest priority
            patient["priority_score"] = 3
            patient["reason"] = "Surgical emergency due to active neurologic deficits (Cauda Equina Syndrome)."
        elif patient["instability"] == "Severe":
            # Severe instability poses a high risk of future neurologic injury
            patient["priority_score"] = 2
            patient["reason"] = "High priority due to severe spinal instability and risk of neurologic injury."
        else:
            # Less severe instability is the lowest priority among these cases
            patient["priority_score"] = 1
            patient["reason"] = "Priority for stabilization, but instability is less severe than other cases."

    # Sort patients by priority score in descending order
    prioritized_list = sorted(patients, key=lambda p: p["priority_score"], reverse=True)

    # Print the final report
    print("Spine Surgery Triage Prioritization Report")
    print("===========================================")
    print("Patients are prioritized from top to lowest priority based on urgency:\n")
    
    for i, patient in enumerate(prioritized_list):
        print(f"Priority Rank {i+1}: Patient {patient['id']}")
        print(f"  Diagnosis: {patient['diagnosis']}")
        print(f"  Reasoning: {patient['reason']}\n")
    
    final_order = [p['id'] for p in prioritized_list]
    print("-------------------------------------------")
    print("The final prioritized order is:")
    print(f"Patient {final_order[0]}, then Patient {final_order[1]}, then Patient {final_order[2]}")
    

solve_triage_task()
<<<F>>>