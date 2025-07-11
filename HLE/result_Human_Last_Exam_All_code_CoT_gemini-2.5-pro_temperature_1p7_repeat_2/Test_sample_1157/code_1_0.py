# Key numerical parameters from the patient's case description.

# Glycated Hemoglobin (HbA1c) indicates long-term blood sugar control.
# A level of 7.5% signifies poorly controlled diabetes, impacting healing.
hba1c = 7.5

# The time delay before receiving hospital care.
# 28 hours is a significant delay, increasing infection risk and ruling out tooth reimplantation.
time_to_treatment_hours = 28

# The SNB angle is a cephalometric measurement in orthodontics and maxillofacial surgery.
# It relates the position of the mandible (lower jaw) to the cranial base.
# A normal angle is ~80 degrees. 88 degrees indicates a prognathic mandible (skeletal Class III).
snb_angle_degrees = 88

# The number of teeth lost in the accident.
teeth_lost = 3
lost_teeth_names = ["Left Lateral Incisor", "Left Canine", "First Premolar"]


# This code formats and prints the summary of the patient's key data.
# This represents a "summary equation" of the factors influencing the treatment plan.
print("--- Patient Case Summary ---")
print(f"Factor 1: Glycated Hemoglobin (HbA1c) = {hba1c} %")
print(f"Factor 2: Time to Treatment = {time_to_treatment_hours} hours")
print(f"Factor 3: SNB Angle = {snb_angle_degrees} degrees")
print(f"Factor 4: Number of Teeth Lost = {teeth_lost}")
print("--------------------------")
