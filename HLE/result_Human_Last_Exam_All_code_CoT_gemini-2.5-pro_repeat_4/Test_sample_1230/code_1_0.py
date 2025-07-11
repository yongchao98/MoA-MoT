import math

# Step 1: Define Charles's known and assumed parameters.
age = 42
bmi = 37
# We don't have his height, so we'll assume an average male height of 175 cm.
assumed_height_cm = 175

# Step 2: Calculate his estimated weight based on BMI and assumed height.
height_m = assumed_height_cm / 100
# Formula: Weight = BMI * (Height in m)^2
estimated_weight_kg = bmi * (height_m ** 2)

print(f"Based on a BMI of {bmi} and an assumed height of {assumed_height_cm} cm, Charles's estimated weight is {estimated_weight_kg:.1f} kg.")

# Step 3: Calculate his Basal Metabolic Rate (BMR) using the Mifflin-St Jeor equation.
# BMR = (10 * weight in kg) + (6.25 * height in cm) - (5 * age) + 5
bmr = (10 * estimated_weight_kg) + (6.25 * assumed_height_cm) - (5 * age) + 5
print(f"His estimated Basal Metabolic Rate (BMR), the calories he burns at rest, is {bmr:.0f} calories per day.")

# Step 4: Calculate his Total Daily Energy Expenditure (TDEE) with a new exercise plan.
# We'll use an activity factor of 1.375 for light exercise (1-3 days/week).
activity_factor = 1.375
tdee = bmr * activity_factor
print(f"With a new exercise plan, his estimated Total Daily Energy Expenditure (TDEE) is {tdee:.0f} calories per day.")

# Step 5: Define a hypothetical intake for his highly restrictive diet and calculate the energy deficit.
# A diet severely limiting protein, fat, AND carbs is a Very Low-Calorie Diet (VLCD).
# We will use a hypothetical value of 1200 calories.
hypothetical_calorie_intake = 1200
energy_deficit = tdee - hypothetical_calorie_intake

print(f"\nA diet restricting all macronutrients is likely very low in calories (e.g., {hypothetical_calorie_intake} calories).")
print("\n--- Calculating the Energy Imbalance ---")
print(f"Required Energy (TDEE): {tdee:.0f} calories")
print(f"Dietary Intake: {hypothetical_calorie_intake} calories")
# Final equation showing each number as requested
print(f"Resulting Daily Deficit = {tdee:.0f} - {hypothetical_calorie_intake} = {energy_deficit:.0f} calories")
print("----------------------------------------")

print("\nConclusion: This large energy deficit is the most likely reason Charles is 'struggling.'")
print("He needs to adjust his diet first to ensure he has enough energy to function and exercise safely.")
print("Specifically, he should not be severely restricting carbohydrates while also limiting fat and protein.")
