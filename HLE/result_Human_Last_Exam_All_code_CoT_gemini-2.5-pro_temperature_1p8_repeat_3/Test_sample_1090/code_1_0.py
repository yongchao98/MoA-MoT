# This script calculates the patient's heart rate response to postural change.

# Patient's heart rate at rest in beats per minute (bpm)
resting_hr = 76

# Patient's heart rate on attempting to stand in bpm
standing_hr = 118

# Calculate the increase in heart rate
hr_increase = standing_hr - resting_hr

print("This patient's primary issue is a failure to progress in physical therapy, making it crucial to identify the limiting factor.")
print("One key objective finding is his heart rate response to changing posture.\n")

print("Calculation of Heart Rate Increase:")
# The final equation prints each number involved in the calculation
print(f"Standing Heart Rate ({standing_hr} bpm) - Resting Heart Rate ({resting_hr} bpm) = Increase of {hr_increase} bpm\n")

print(f"An increase of {hr_increase} bpm is significant and indicates postural tachycardia. This, combined with the patient's low BMI of 18.5 kg/m2, points to severe deconditioning and malnutrition as the most likely reasons for his inability to ambulate.")
print("Therefore, the most appropriate next step is to address his nutritional status to provide the energy and protein required for recovery and successful rehabilitation.")