# Patient's vital signs
hr_rest = 76  # Heart rate at rest (beats per minute)
hr_stand = 118 # Heart rate when standing (beats per minute)

# The threshold for postural tachycardia is a heart rate increase of 30 bpm or more upon standing.
pots_threshold = 30

# Calculate the change in heart rate
hr_change = hr_stand - hr_rest

# Print the analysis
print("The patient's clinical picture suggests severe deconditioning and malnutrition are impeding his recovery.")
print("One key finding is his heart rate response to standing.")
print(f"The equation for the change in heart rate is: {hr_stand} (standing) - {hr_rest} (resting) = {hr_change} bpm.")

if hr_change >= pots_threshold:
    print(f"The heart rate increases by {hr_change} bpm, which meets the criteria for postural tachycardia (>= {pots_threshold} bpm).")
    print("This finding often points to deconditioning and volume depletion, which can be addressed with nutritional and fluid support.")
else:
    print(f"The heart rate increase of {hr_change} bpm does not meet the criteria for postural tachycardia.")

# Note on BMI
bmi = 18.5
print(f"\nAdditionally, his BMI is {bmi} kg/m2, which is borderline underweight and a major risk factor for poor outcomes and failure to recover from illness.")