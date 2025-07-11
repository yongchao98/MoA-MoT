def provide_diagnosis():
    """
    This function outlines the reasoning for the medical diagnosis based on the provided case study.
    """
    patient_profile = {
        "Age": "62-year-old male",
        "History": "20-pack-year smoker, shipyard worker"
    }
    
    initial_symptoms = [
        "Fatigue",
        "Polyarthritis (wrists, ankles, elbows)",
        "Multiple pulmonary nodules"
    ]
    
    progressive_symptoms = [
        "Dizziness and confusion (neurological)",
        "Bruising (vasculitic sign)",
        "Difficulty swallowing",
        "Shortness of breath"
    ]
    
    terminal_event = [
        "Immunosuppression from steroid therapy",
        "Secondary infection (fever, productive cough)",
        "Death from septic shock"
    ]
    
    primary_diagnosis = "Granulomatosis with Polyangiitis (GPA)"
    
    print("Based on the analysis of the case, the most likely underlying disease is:")
    print(f"*** {primary_diagnosis} ***")
    print("\nReasoning:")
    print("1. The combination of lung involvement (multiple nodules) and systemic symptoms (polyarthritis) points to a systemic disease.")
    print("2. The progressive symptoms, including neurological changes and bruising, are characteristic of a small-to-medium vessel vasculitis, a hallmark of GPA.")
    print("3. The patient's death from septic shock is a known, serious complication for individuals with GPA who are on immunosuppressive therapy (like steroids), making them vulnerable to severe infections.")

provide_diagnosis()