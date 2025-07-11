def analyze_lifestyle_change(age, bmi, diet_details, exercise_details):
    """
    Analyzes a new lifestyle plan to suggest the first point of adjustment.
    """
    # Assign difficulty scores to represent the intensity and restrictiveness of the changes.
    # A diet limiting three major macronutrients is exceptionally difficult.
    diet_difficulty = 9  # on a scale of 1-10
    # Starting a new exercise routine while traveling is also challenging.
    exercise_difficulty = 7 # on a scale of 1-10

    # The problem states Charles is "struggling".
    is_struggling = True
    
    # The total demand from making two difficult changes at once is very high.
    total_demand = diet_difficulty + exercise_difficulty

    first_adjustment = "None"
    reasoning = ""

    # A sustainable plan is crucial. If the total demand is too high and causing struggles,
    # the most impactful and difficult component should be adjusted first.
    if is_struggling and total_demand > 10: # A subjective threshold for "too much at once"
        if diet_difficulty > exercise_difficulty:
            first_adjustment = "the diet"
            reasoning = (
                "The diet is extremely restrictive. Limiting carbohydrates, fats, and animal protein "
                "simultaneously can severely deplete energy levels, making it nearly impossible to "
                "start and maintain a new exercise program."
            )
        else:
            first_adjustment = "the exercise"
            reasoning = (
                "The exercise routine is too intense to start with. It's better to build a solid "
                "dietary foundation first before adding strenuous physical activity."
            )
    
    print(f"Analysis for Charles, a {age}-year-old man with a BMI of {bmi}:")
    print(f"He has started a demanding new plan and is struggling.")
    print(f"His diet is '{diet_details}' (Difficulty Score: {diet_difficulty}).")
    print(f"His exercise plan is '{exercise_details}' (Difficulty Score: {exercise_difficulty}).")
    print("-" * 20)
    print(f"Recommendation: Charles needs to adjust {first_adjustment} first.")
    print(f"\nJustification: {reasoning}")
    print("He should aim for a more balanced intake, particularly ensuring enough complex carbohydrates for energy, before increasing his exercise intensity.")


# Run the analysis for Charles's situation.
charles_age = 42
charles_bmi = 37
charles_diet = "Whole foods, limiting animal protein, fat, and carbohydrates"
charles_exercise = "New workout routine in hotel gyms"

analyze_lifestyle_change(charles_age, charles_bmi, charles_diet, charles_exercise)