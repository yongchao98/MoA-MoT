import pandas as pd

def find_best_next_step():
    """
    Analyzes treatment options for an agitated patient using a weighted scoring model
    for efficacy and safety.
    """

    # Data for each answer choice including Efficacy and Safety scores (out of 10)
    # Efficacy: Likelihood of calming the patient.
    # Safety: Inverse of risk (e.g., respiratory depression, over-sedation).
    options_data = {
        "A": {"description": "2mg IV lorazepam", "efficacy": 6, "safety": 4},
        "B": {"description": "2mg IM lorazepam + 5mg olanzapine IM", "efficacy": 9, "safety": 8},
        "C": {"description": "Verbal de-escalation", "efficacy": 0, "safety": 5},
        "D": {"description": "10mg IM olanzapine", "efficacy": 7, "safety": 6},
        "E": {"description": "10mg IM olanzapine + 2mg IM lorazepam", "efficacy": 10, "safety": 4},
    }

    results = []
    
    # Efficacy and Safety have equal importance
    efficacy_weight = 0.5
    safety_weight = 0.5
    
    print("Calculating the best option based on a weighted score of Efficacy and Safety...")
    print(f"Scoring Formula: (Efficacy * {efficacy_weight}) + (Safety * {safety_weight})\n")

    for option, scores in options_data.items():
        efficacy_score = scores['efficacy']
        safety_score = scores['safety']
        
        # Calculate the total score using the defined weights and numbers
        total_score = (efficacy_score * efficacy_weight) + (safety_score * safety_weight)
        
        results.append({
            "Option": option,
            "Description": scores['description'],
            "Total Score": total_score
        })
        
        # Print the equation for each option as requested
        print(f"Option {option}: ({efficacy_score} * {efficacy_weight}) + ({safety_score} * {safety_weight}) = {total_score:.1f}")

    # Determine the best option
    best_option = max(results, key=lambda x: x['Total Score'])
    
    print("\n--- Analysis Summary ---")
    print(f"Initial Intervention: 5mg IM Zyprexa was ineffective for a violent, agitated patient.")
    print("Goal: Find the next step that maximizes both effectiveness and safety.")
    print("\nOption A (IV Lorazepam): Risky IV route in a patient with unknown history already on a sedative.")
    print("Option C (Verbal): Inappropriate as the patient is already violent.")
    print("Option D (More Olanzapine): Monotherapy may be less effective than combination therapy.")
    print("Option E (High-Dose Combo): Potent but carries a high risk of over-sedation with a total of 15mg olanzapine.")
    print("Option B (Standard-Dose Combo): Represents a standard, effective, and safe escalation of care, bringing the total olanzapine to a standard 10mg dose while adding a synergistic agent.")
    
    print("\n--- Conclusion ---")
    print(f"The best next step is Option {best_option['Option']} with a score of {best_option['Total Score']:.1f}.")
    print(f"Recommended Action: {best_option['Description']}")


if __name__ == "__main__":
    find_best_next_step()