# Patient's relevant clinical data
bmi = 18.5  # in kg/m2
resting_hr = 76  # in beats/min
standing_hr = 118  # in beats/min

# --- Analysis ---

# 1. Analyze Body Mass Index (BMI)
# A BMI of 18.5 is on the threshold between normal and underweight.
# In an elderly, hospitalized patient, this strongly suggests a poor nutritional state and sarcopenia (muscle loss).
print(f"Patient's BMI is {bmi} kg/m2, which indicates a borderline underweight/poor nutritional status.")

# 2. Calculate and analyze the heart rate response to postural change
hr_increase = standing_hr - resting_hr
print(f"\nThe heart rate increases from {resting_hr}/min at rest to {standing_hr}/min upon standing.")
print(f"This is an increase of {standing_hr} - {resting_hr} = {hr_increase} beats/min.")

# 3. Interpretation
print("\n--- Interpretation ---")
print("A heart rate increase of over 30 beats/min with minimal exertion like standing is a sign of significant cardiovascular deconditioning.")
print("This profound deconditioning, combined with the low BMI, is the most likely reason the patient cannot participate effectively in physical therapy.")
print("To rebuild muscle and endurance, the patient requires adequate caloric and protein intake.")
print("\nTherefore, addressing his nutritional deficit is the most critical next step to enable recovery and successful rehabilitation.")
