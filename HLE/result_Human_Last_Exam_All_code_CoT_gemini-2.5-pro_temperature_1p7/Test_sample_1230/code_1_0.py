import math

def analyze_charles_diet():
    """
    Analyzes Charles's dietary situation to identify the primary issue.
    """
    # Step 1: Define knowns and make reasonable assumptions
    age = 42
    bmi = 37
    # Assume a standard height of 175cm (approx. 5'9") to calculate weight
    height_cm = 175
    height_m = height_cm / 100

    # Step 2: Calculate estimated weight from BMI
    # Formula: weight(kg) = BMI * height(m)^2
    weight_kg = bmi * (height_m ** 2)

    # Step 3: Calculate Basal Metabolic Rate (BMR) using Mifflin-St Jeorge formula
    # Formula: BMR = 10 * weight(kg) + 6.25 * height(cm) - 5 * age + 5 (for men)
    bmr = (10 * weight_kg) + (6.25 * height_cm) - (5 * age) + 5

    # Step 4: Estimate Total Daily Energy Expenditure (TDEE)
    # Assume a "lightly active" factor as he plans to start exercise
    activity_factor = 1.375
    tdee = bmr * activity_factor

    # Step 5: Model the overly restrictive diet and compare
    # A diet limiting all macronutrients is very low calorie. Assume 1200 kcal.
    struggling_diet_kcal = 1200
    current_deficit = tdee - struggling_diet_kcal
    
    # A sustainable deficit is around 500 kcal for 1 lb/week weight loss
    recommended_deficit = 500
    recommended_diet_kcal = tdee - recommended_deficit

    # Step 6: Print the analysis
    print("--- Charles's Energy Needs Analysis ---")
    print("\n[1] Estimating Physical Statistics:")
    print(f"Based on Age={age}, BMI={bmi}, and an assumed height of {height_cm}cm...")
    print(f"Estimated Weight = {bmi} * {height_m:.2f} * {height_m:.2f} = {weight_kg:.1f} kg")

    print("\n[2] Calculating Daily Calorie Needs:")
    print("Using the Mifflin-St Jeorge formula for Basal Metabolic Rate (BMR):")
    print(f"BMR = (10 * {weight_kg:.1f}) + (6.25 * {height_cm}) - (5 * {age}) + 5 = {math.ceil(bmr)} calories/day")
    print("\nTo maintain his current weight with light activity, his Total Daily Energy Expenditure (TDEE) is:")
    print(f"TDEE = {math.ceil(bmr)} (BMR) * {activity_factor} (Activity Factor) = {math.ceil(tdee)} calories/day")
    
    print("\n[3] Analyzing the Current Diet's Impact:")
    print(f"An extremely restrictive diet might provide only ~{struggling_diet_kcal} calories.")
    print("This creates a massive daily energy deficit:")
    print(f"Deficit = {math.ceil(tdee)} (TDEE) - {struggling_diet_kcal} (Intake) = {math.ceil(current_deficit)} calories")
    
    print("\n--- Conclusion ---")
    print(f"A sustainable deficit for healthy weight loss is typically around {recommended_deficit} calories.")
    print(f"Charles's current deficit of {math.ceil(current_deficit)} is too large, leading to fatigue and making the diet unsustainable.")
    print("The first and most important adjustment Charles needs to make is to his diet.")
    print(f"He should adopt a more balanced eating plan with a moderate deficit, aiming for around {math.ceil(recommended_diet_kcal)} calories, not severely limiting all macronutrients.")

analyze_charles_diet()