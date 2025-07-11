def analyze_charles_diet_plan():
    """
    Calculates Charles's estimated caloric needs to show why his diet is too restrictive.
    """
    # Charles's known information
    age = 42
    bmi = 37

    # Step 1: Assume a reasonable height to estimate his weight.
    # We'll assume an average height for a man: 1.78 meters (approx. 5'10").
    height_m = 1.78
    height_cm = height_m * 100

    # Calculate weight from BMI (BMI = weight_kg / height_m^2)
    weight_kg = bmi * (height_m ** 2)

    print(f"Charles's Profile:")
    print(f"- Age: {age} years")
    print(f"- BMI: {bmi}")
    print(f"- Assumed Height: {height_m} m")
    print(f"- Estimated Weight: {weight_kg:.1f} kg\n")

    # Step 2: Calculate Basal Metabolic Rate (BMR) using the Mifflin-St Jeor equation.
    # BMR = 10 * weight (kg) + 6.25 * height (cm) - 5 * age (years) + 5 (for men)
    bmr = (10 * weight_kg) + (6.25 * height_cm) - (5 * age) + 5
    
    print("Calculating Basal Metabolic Rate (BMR):")
    # Outputting each number in the final equation as requested
    print(f"BMR = (10 * {weight_kg:.1f}) + (6.25 * {height_cm:.0f}) - (5 * {age}) + 5")
    print(f"BMR = {10 * weight_kg:.1f} + {6.25 * height_cm:.1f} - {5 * age} + 5 = {bmr:.1f} calories/day\n")
    print("This is the number of calories his body burns at rest.\n")

    # Step 3: Calculate Total Daily Energy Expenditure (TDEE).
    # Since he is starting to exercise, we'll use a 'lightly active' multiplier.
    activity_multiplier = 1.375 # Light exercise (1-3 days/week)
    tdee = bmr * activity_multiplier

    print("Calculating Total Daily Energy Expenditure (TDEE):")
    print(f"TDEE = BMR * Activity Multiplier")
    print(f"TDEE = {bmr:.1f} * {activity_multiplier} = {tdee:.1f} calories/day\n")
    print(f"This is the estimated number of calories Charles needs to maintain his current weight while lightly exercising.\n")

    # Step 4: Suggest a sustainable calorie target for weight loss.
    # A 500-750 calorie deficit per day is a standard, sustainable goal.
    calorie_deficit = 750
    weight_loss_target = tdee - calorie_deficit
    
    print("Conclusion & Recommendation:")
    print(f"To lose weight sustainably, a common goal is a {calorie_deficit} calorie deficit.")
    print(f"Recommended Target: {tdee:.1f} - {calorie_deficit} = {weight_loss_target:.1f} calories/day.")
    print("\nCharles is struggling because his diet, which limits protein, fat, AND carbs, is likely providing far fewer calories than this target.")
    print("This extreme restriction leads to low energy, making both the diet and exercise feel impossible.")
    print("\nTherefore, the first thing Charles needs to adjust is his overly restrictive diet to ensure he has enough energy to function and exercise.")

analyze_charles_diet_plan()