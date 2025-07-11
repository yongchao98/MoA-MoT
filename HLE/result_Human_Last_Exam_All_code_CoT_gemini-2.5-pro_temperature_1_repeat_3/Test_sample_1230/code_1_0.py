def create_plan_adjustment_summary():
    """
    Analyzes Charles's situation and identifies the primary area for adjustment.
    """
    name = "Charles"
    age = 42
    bmi = 37
    primary_issue = "Struggling to keep up with self-imposed demands."
    diet_details = "Limiting animal protein, fat, AND carbohydrates simultaneously."
    exercise_details = "Started a new exercise plan."
    lifestyle = "Travels frequently for work."

    # Analysis
    analysis = (
        f"{name}, a {age}-year-old with a BMI of {bmi}, is finding his new plan challenging. "
        f"The core problem lies in the extreme nature of his diet. A diet that severely restricts all three macronutrients "
        f"(protein, fat, and carbohydrates) is not sustainable and leads to very low energy levels. "
        f"This lack of fuel makes it nearly impossible to begin a new exercise routine and manage a demanding travel schedule."
    )

    # Recommendation
    recommendation = (
        "The first and most critical adjustment Charles needs to make is to his diet. "
        "He should move away from restricting everything at once and adopt a more balanced approach. "
        "He needs to ensure he is consuming enough complex carbohydrates for energy, adequate protein for muscle repair, and healthy fats for overall function. "
        "Adjusting the diet to be more sustainable will provide the foundation he needs to successfully incorporate exercise and achieve his long-term health goals."
    )

    print("--- Analysis of Charles's Situation ---")
    print(analysis)
    print("\n--- Recommended First Adjustment ---")
    print(recommendation)

create_plan_adjustment_summary()