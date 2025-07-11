def analyze_health_plan():
    """
    Analyzes Charles's health plan to identify the most critical area for adjustment.
    """
    # Charles's Profile Data
    age = 42
    bmi = 37

    # Plan Components Analysis
    # The diet plan is the most significant immediate change.
    diet_limits_animal_protein = True
    diet_limits_fat = True
    diet_limits_carbohydrates = True
    
    # Psychotherapy is reported as helpful.
    psychotherapy_is_helpful = True
    
    # Exercise is still in the planning phase.
    exercise_started = False
    
    print("Analyzing Charles's Health Plan:")
    print(f"Age: {age}")
    print(f"BMI: {bmi}")
    print("-" * 30)

    # The core logic: check for an overly restrictive diet.
    # Limiting all three macronutrients is a common cause of failure in new diets.
    if diet_limits_fat and diet_limits_carbohydrates and diet_limits_animal_protein:
        print("Diagnosis: The current diet is extremely restrictive.")
        print("Reasoning: Limiting carbohydrates, fats, and protein sources all at once can lead to severe energy deficit, fatigue, and hunger.")
        print("This is the most likely reason Charles is 'struggling'.")
        print("\nRecommendation: The first adjustment Charles needs to make is to his diet.")
        adjustment_needed = "His diet is too restrictive and needs to be more balanced."
    elif not psychotherapy_is_helpful:
        adjustment_needed = "His psychotherapy approach may need to be re-evaluated."
    else:
        adjustment_needed = "The plan seems balanced, more information is needed on the struggle."

    print("\n--- Conclusion ---")
    print(f"Where does Charles need to adjust first? \n> {adjustment_needed}")


analyze_health_plan()