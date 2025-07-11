# The patient's clinical vignette provides key vital signs that help identify the primary problem.
# Let's analyze the change in heart rate from resting to standing.

# Resting heart rate (in beats per minute)
hr_rest = 76

# Heart rate on attempting to stand (in beats per minute)
hr_standing = 118

# Calculate the increase in heart rate
hr_increase = hr_standing - hr_rest

# This calculation quantifies the patient's significant postural tachycardia,
# which is a sign of severe deconditioning. Addressing this functional
# impairment is the most critical next step.
print(f"The patient's heart rate increased by {hr_standing} - {hr_rest} = {hr_increase} beats per minute upon standing.")