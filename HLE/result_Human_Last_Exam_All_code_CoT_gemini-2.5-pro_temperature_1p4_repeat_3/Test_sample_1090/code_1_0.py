# Patient's heart rate data
resting_hr = 76  # beats/min
standing_hr = 118 # beats/min

# Calculate the change in heart rate
hr_change = standing_hr - resting_hr

# Print the analysis
print("Analyzing the patient's heart rate response to postural change:")
print(f"Resting Heart Rate: {resting_hr} beats/min")
print(f"Heart Rate on Standing: {standing_hr} beats/min")
print("The change in heart rate is calculated as the standing heart rate minus the resting heart rate.")
print(f"Calculation: {standing_hr} - {resting_hr} = {hr_change} beats/min")
print(f"\nA heart rate increase of {hr_change} beats/min upon standing is significant and indicates severe deconditioning, which explains the patient's inability to tolerate activity.")
print("This supports the need for intensive, specialized rehabilitation.")