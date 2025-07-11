# The patient has a new, significant functional impairment that explains his
# inability to walk and points towards the necessary next step in his care.
# This script calculates the severity of the knee flexion contracture, which
# is a major mechanical barrier to ambulation.

full_knee_extension_degrees = 180
patient_knee_extension_degrees = 165

# A flexion contracture is the inability to fully straighten or extend a joint.
# The severity is the difference between full extension (180 degrees for the knee)
# and the patient's maximum extension.
flexion_contracture_degrees = full_knee_extension_degrees - patient_knee_extension_degrees

print("The patient is medically stabilizing but cannot walk due to severe deconditioning and a new physical limitation.")
print(f"A key finding is the inability to fully extend the left knee to the normal {full_knee_extension_degrees} degrees.")
print(f"The patient's maximum knee extension is only {patient_knee_extension_degrees} degrees.")
print(f"This represents a flexion contracture of: {full_knee_extension_degrees} - {patient_knee_extension_degrees} = {flexion_contracture_degrees} degrees.")
print("\nThis new contracture, combined with his deconditioning and poor nutritional status (BMI 18.5), requires intensive, multidisciplinary therapy that is best provided in a specialized setting.")
print("\nTherefore, the single most appropriate next step in management is to arrange this specialized care.")
