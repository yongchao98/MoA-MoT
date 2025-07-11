import pandas as pd

def prioritize_patients():
    """
    Prioritizes spine surgery patients based on clinical presentation.
    The function uses a scoring system where neurological deficits are weighted
    most heavily, followed by the severity of the fracture.
    """
    
    # Define patient data
    patients = [
        {"id": 1, "condition": "Severe burst fracture", "neuro_deficit": False, "neuro_description": "no neurologic deficits"},
        {"id": 2, "condition": "Compression fracture with mild spondylolisthesis", "neuro_deficit": False, "neuro_description": "no neurologic deficits"},
        {"id": 3, "condition": "Split fracture", "neuro_deficit": True, "neuro_description": "mildly disordered pelvic functions"}
    ]

    # Define scoring criteria
    # Neurological deficit is a critical factor and receives a very high score.
    neuro_deficit_score = 100
    
    # Fracture severity contributes to the score based on instability.
    fracture_scores = {
        "Severe burst fracture": 50,
        "Split fracture": 30,
        "Compression fracture with mild spondylolisthesis": 10
    }

    print("Prioritizing patients based on a clinical scoring system.\n")
    print("Priority Score = (Neurological Deficit Score) + (Fracture Severity Score)\n")
    
    scored_patients = []
    for patient in patients:
        # Calculate the score from the neurological status
        neuro_score = neuro_deficit_score if patient["neuro_deficit"] else 0
        
        # Calculate the score from the fracture type
        fracture_score = fracture_scores.get(patient["condition"], 0)
        
        # Calculate the total priority score
        total_score = neuro_score + fracture_score
        
        scored_patients.append({
            "id": patient["id"],
            "total_score": total_score,
            "neuro_score": neuro_score,
            "fracture_score": fracture_score,
            "neuro_description": patient["neuro_description"]
        })

        # Output the "equation" for each patient's score
        print(f"Patient {patient['id']} Priority = {neuro_score} (due to: {patient['neuro_description']}) + {fracture_score} (due to: {patient['condition']}) = {total_score}")


    # Sort patients by score in descending order (highest score is highest priority)
    sorted_patients = sorted(scored_patients, key=lambda x: x["total_score"], reverse=True)
    
    print("\n--- Final Prioritized List (from highest to lowest priority) ---")
    priority_order = []
    for i, patient in enumerate(sorted_patients):
        priority_order.append(str(patient['id']))
        print(f"{i+1}. Patient {patient['id']} (Score: {patient['total_score']})")

    final_order_string = ", ".join(priority_order)
    print(f"\nThe final priority order is: Patient {final_order_string}.")

prioritize_patients()
<<<F>>>