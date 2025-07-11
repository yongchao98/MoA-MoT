# This script converts the newborn's weight from imperial to metric units.

# Given weight
pounds = 12
ounces = 1

# Conversion factors
# 1 pound = 0.453592 kilograms
# 1 ounce = 0.0283495 kilograms
lb_to_kg = 0.453592
oz_to_kg = 0.0283495

# Calculate total weight in kilograms
total_kg = (pounds * lb_to_kg) + (ounces * oz_to_kg)

# The final equation is: (12 * 0.453592) + (1 * 0.0283495) = total_kg
# The numbers in the final equation are 12, 1, and the result.
# We will print the input values and the final result.

print(f"Patient's weight in pounds: {pounds}")
print(f"Patient's weight in ounces: {ounces}")
print(f"The newborn's total weight is approximately {total_kg:.2f} kg, which is classified as macrosomia (large for gestational age).")
