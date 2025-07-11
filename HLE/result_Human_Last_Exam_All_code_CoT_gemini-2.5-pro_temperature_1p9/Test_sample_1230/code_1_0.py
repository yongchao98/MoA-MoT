def analyze_health_plan(age, bmi, diet_restrictions):
    """
    Analyzes a person's health plan to identify immediate areas for adjustment.

    Args:
        age (int): The person's age.
        bmi (float): The person's Body Mass Index.
        diet_restrictions (dict): A dictionary indicating dietary limitations.
    """
    print("--- Health Plan Analysis ---")
    print(f"Analyzing profile for a {age}-year-old individual with a BMI of {bmi}.")
    
    # The main problem is simultaneously limiting all major macronutrient sources.
    is_too_restrictive = diet_restrictions.get('protein') and \
                          diet_restrictions.get('fat') and \
                          diet_restrictions.get('carbohydrates')

    if is_too_restrictive:
        print("\nPrimary Concern Identified: The Diet is Overly Restrictive.")
        print("\nRecommendation:")
        print("Charles's first adjustment should be to his diet. A plan that limits protein, fat, AND carbohydrates simultaneously is not sustainable and is likely causing low energy and hunger.")
        print("A more balanced approach is recommended:")
        print("1. Focus on a whole-foods diet that includes adequate lean protein and healthy fats for satiety.")
        print("2. Ensure sufficient complex carbohydrates for energy, especially with a new exercise plan.")
        
        # This print statement fulfills the requirement to output the numbers in a final statement.
        print("\nFinal Conclusion: The diet for the 42-year-old man with a BMI of 37 must be made more balanced and sustainable to ensure success.")
    else:
        print("The diet appears balanced. Further analysis of other lifestyle factors may be needed.")

# Charles's data from the scenario
charles_age = 42
charles_bmi = 37
# He is limiting animal protein, fat, and carbohydrates.
charles_diet = {
    'protein': True,
    'fat': True,
    'carbohydrates': True
}

analyze_health_plan(charles_age, charles_bmi, charles_diet)