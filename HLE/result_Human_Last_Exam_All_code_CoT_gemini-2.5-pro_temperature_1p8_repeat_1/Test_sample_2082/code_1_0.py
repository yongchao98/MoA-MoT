# The patient's data from the case study
weight_kg = 65.77
height_m = 1.68

# Calculate the Body Mass Index (BMI)
# The formula is weight (kg) / (height (m))^2
bmi = weight_kg / (height_m ** 2)

# Define BMI categories
if bmi < 18.5:
    category = "Underweight"
elif 18.5 <= bmi < 24.9:
    category = "Normal weight"
elif 25 <= bmi < 29.9:
    category = "Overweight"
else:
    category = "Obesity"

# Output the calculation and the result
print("Calculating the patient's Body Mass Index (BMI):")
print(f"BMI = weight / (height * height)")
print(f"BMI = {weight_kg} kg / ({height_m} m * {height_m} m)")
print(f"The calculated BMI is: {bmi:.2f}")
print(f"This BMI falls into the '{category}' category.")