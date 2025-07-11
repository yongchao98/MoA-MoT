def solve_charles_diet_problem():
    """
    Analyzes Charles's dietary problem by calculating his energy needs
    and showing the unsustainability of his current diet plan.
    """
    # Given data
    age = 42  # years
    bmi = 37

    # Step 1: Estimate Weight from BMI. Assume an average height for a man.
    # We'll assume a height of 1.80 meters (approx. 5'11").
    height_m = 1.80
    height_cm = height_m * 100
    # BMI = weight (kg) / height (m)^2  => weight = BMI * height^2
    weight_kg = bmi * (height_m ** 2)

    print(f"Assuming a height of {height_m}m for a BMI of {bmi}, estimated weight is {weight_kg:.1f}kg.")

    # Step 2: Calculate Basal Metabolic Rate (BMR) using Mifflin-St Jeor equation for men.
    # BMR = 10 * weight (kg) + 6.25 * height (cm) - 5 * age (years) + 5
    bmr = (10 * weight_kg) + (6.25 * height_cm) - (5 * age) + 5
    print(f"His estimated Basal Metabolic Rate (BMR) is {bmr:.0f} kcal/day.")

    # Step 3: Calculate Total Daily Energy Expenditure (TDEE).
    # Since he started exercising, we'll use a "lightly active" multiplier.
    activity_factor = 1.375
    tdee = bmr * activity_factor
    print(f"His estimated Total Daily Energy Expenditure (TDEE) with light exercise is {tdee:.0f} kcal/day.")

    # Step 4: Calculate a target calorie intake for weight loss (a 500 kcal deficit is common).
    calorie_deficit = 500
    target_calories = tdee - calorie_deficit
    print(f"A sustainable weight loss plan could target {target_calories:.0f} kcal/day.")
    print("-" * 30)

    # Step 5: Model his overly restrictive diet (low-carb AND low-fat).
    # This demonstrates why he is struggling.
    # Assume "limited" means very low intake for both.
    restricted_carbs_g = 50
    restricted_fat_g = 30

    calories_from_carbs = restricted_carbs_g * 4
    calories_from_fat = restricted_fat_g * 9

    print("Analyzing the overly restrictive diet (limiting both carbs and fat):")
    print(f"Calories from {restricted_carbs_g}g of Carbs = {calories_from_carbs} kcal")
    print(f"Calories from {restricted_fat_g}g of Fat = {calories_from_fat} kcal")
    
    # Step 6: Calculate the required protein to meet the target, which will be unsustainable.
    required_protein_calories = target_calories - calories_from_carbs - calories_from_fat
    required_protein_grams = required_protein_calories / 4
    
    print("\nTo meet his calorie target, his protein requirement would be:")
    # The final equation requested by the user prompt
    print(f"Required Protein Grams = ({target_calories:.0f} - {calories_from_carbs} - {calories_from_fat}) / 4")
    print(f"Required Protein Grams = {required_protein_grams:.1f}g")

    print("\nConclusion: Consuming over {required_protein_grams:.0f}g of protein, while also limiting animal protein, is unsustainable.".format(required_protein_grams=required_protein_grams))
    print("The primary issue is the diet's composition. It's too restrictive on all macronutrients.")
    print("Charles needs to adjust his diet to ensure he has enough energy, for instance, by moderately increasing either healthy fats or complex carbohydrates.")

solve_charles_diet_problem()