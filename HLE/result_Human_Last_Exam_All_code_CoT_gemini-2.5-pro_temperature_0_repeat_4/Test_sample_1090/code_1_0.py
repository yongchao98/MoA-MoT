# The patient's Body Mass Index (BMI) is a key indicator of his poor nutritional status,
# which is a major contributor to his weakness and inability to participate in rehabilitation.
# BMI is calculated as weight in kilograms divided by the square of height in meters.
# The problem states the BMI is 18.5 kg/m2. Let's demonstrate the calculation.
# We can assume a height for the patient, for example, 1.78 meters (approx. 5'10").

# Given values
bmi = 18.5  # kg/m^2
height_m = 1.78  # meters

# To show the calculation, we first need to find the corresponding weight.
# Weight = BMI * (Height^2)
weight_kg = bmi * (height_m ** 2)

# Now we can demonstrate the BMI calculation with the derived weight and assumed height.
print("Calculating Body Mass Index (BMI):")
print(f"Patient's assumed height: {height_m} m")
print(f"Patient's calculated weight: {round(weight_kg, 1)} kg")
print("-" * 30)
print("BMI Formula: Weight (kg) / [Height (m)]^2")
print(f"Calculation: {round(weight_kg, 1)} / ({height_m} * {height_m})")

# Perform the final calculation
calculated_bmi = round(weight_kg, 1) / (height_m * height_m)

print(f"Result: {round(calculated_bmi, 1)} kg/m^2")
print("\nThis low BMI highlights the patient's poor nutritional status, a key factor that must be addressed.")
