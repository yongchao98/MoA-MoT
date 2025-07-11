# This script calculates the patient's heart rate change on standing.
# This significant change highlights his level of deconditioning.

# Patient's heart rate data from the clinical case
heart_rate_at_rest = 76 # beats per minute
heart_rate_on_standing = 118 # beats per minute

# Calculate the increase in heart rate
heart_rate_increase = heart_rate_on_standing - heart_rate_at_rest

# Print the final equation showing the calculation
print("The patient's postural heart rate change calculation:")
print(f"{heart_rate_on_standing} bpm (standing) - {heart_rate_at_rest} bpm (resting) = {heart_rate_increase} bpm increase")
