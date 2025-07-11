# Plan:
# 1. Define the patient's resting and standing heart rates from the medical case.
# 2. Calculate the increase in heart rate when the patient attempts to stand.
# 3. Print the initial values, the calculated increase, and a clinical interpretation.
# This calculation will quantify the patient's severe deconditioning, which points to a foundational issue like malnutrition hindering recovery.

# Patient's heart rate data
resting_hr = 76 # beats/min
standing_hr = 118 # beats/min

# Calculate the change in heart rate
hr_increase = standing_hr - resting_hr

# Output the results and interpretation
print(f"Patient's resting heart rate: {resting_hr}/min")
print(f"Patient's heart rate on attempting to stand: {standing_hr}/min")
print(f"The increase in heart rate is {standing_hr} - {resting_hr} = {hr_increase} beats/min.")
print("\nThis large increase in heart rate with a simple postural change indicates significant physiological stress and cardiovascular deconditioning.")
print("This condition is severely exacerbated by the patient's catabolic state and borderline low BMI of 18.5 kg/m2.")
print("Therefore, providing the necessary fuel for the body is the most critical next step to enable recovery and tolerate physical therapy.")
