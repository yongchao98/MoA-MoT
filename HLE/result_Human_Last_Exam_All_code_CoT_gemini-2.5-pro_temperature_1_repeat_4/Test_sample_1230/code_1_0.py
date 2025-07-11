import math

# Charles's known data
age = 42  # years
bmi = 37

# Assumed data for calculation
height_m = 1.80  # meters (assuming a height of ~5'11")

# --- Step 1: Calculate weight from BMI ---
# Formula: weight (kg) = BMI * height (m)^2
weight_kg = bmi * (height_m ** 2)
height_cm = height_m * 100

print(f"Charles's Data:")
print(f"Age: {age} years")
print(f"BMI: {bmi}")
print(f"Assumed Height: {height_m:.2f} m")
print(f"Calculated Weight: {weight_kg:.2f} kg\n")


# --- Step 2: Calculate Basal Metabolic Rate (BMR) using Mifflin-St Jeor equation ---
# BMR = (10 * weight in kg) + (6.25 * height in cm) - (5 * age) + 5
print("Calculating Basal Metabolic Rate (BMR) to estimate daily calorie needs at rest:")
print("BMR = (10 * weight) + (6.25 * height) - (5 * age) + 5")

# Breaking down the equation to show each number
weight_term = 10 * weight_kg
height_term = 6.25 * height_cm
age_term = 5 * age
male_adjustment = 5

bmr = weight_term + height_term - age_term + male_adjustment

print(f"The parts of the equation are:")
print(f"  10 * {weight_kg:.2f} = {weight_term:.2f}")
print(f"  6.25 * {height_cm:.2f} = {height_term:.2f}")
print(f"  5 * {age} = {age_term:.2f}")
print(f"  Male adjustment = {male_adjustment}")
print(f"Final Equation: {weight_term:.2f} + {height_term:.2f} - {age_term:.2f} + {male_adjustment} = {bmr:.2f}\n")


# --- Step 3: Analysis and Conclusion ---
print("--- Analysis ---")
print(f"Charles's estimated BMR is {math.ceil(bmr)} calories per day.")
print("This is the energy his body needs just to function at rest.")
print("His total daily energy needs, including his travel and new exercise routine, are significantly higher.")

print("\nA diet that severely limits protein, fat, AND carbohydrates is extremely low in calories,")
print("likely providing far less energy than his BMR. This massive energy deficit would lead to fatigue,")
print("hunger, and difficulty sustaining daily activities, let alone a new workout regimen.\n")

print("--- Conclusion ---")
print("Charles is struggling because his diet is too restrictive and not providing enough fuel for his body.")
print("Therefore, the first and most critical area he needs to adjust is his diet.")
print("He should aim for a more balanced approach that creates a moderate calorie deficit without")
print("drastically eliminating essential macronutrients like protein and healthy fats.")

# To provide a clear answer based on potential choices:
# A. His exercise routine
# B. His travel schedule
# C. His diet
# D. His psychotherapy sessions
# The clear answer is C.
