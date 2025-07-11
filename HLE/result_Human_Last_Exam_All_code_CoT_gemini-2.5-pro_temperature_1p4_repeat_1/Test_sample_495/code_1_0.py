import math

# Given values
E0 = 8.5  # Initial energy in MeV
R0 = 8.3  # Range in cm
x = 4.0   # Distance from the source in cm

# Calculate the energy loss per centimeter using the derived formula
# Formula: |dE/dx| = (2 * E0) / (3 * R0 * (1 - x / R0)^(1/3))

# Calculate the denominator term
denominator_term = (1 - x / R0)**(1/3)

# Calculate the final energy loss per cm
dEdx = (2 * E0) / (3 * R0 * denominator_term)

# Print the final result, showing the equation with the numbers plugged in
print("The formula for the energy loss per centimeter |dE/dx| at a distance x is:")
print("|dE/dx| = (2 * E0) / (3 * R0 * (1 - x / R0)^(1/3))")
print("\nPlugging in the values:")
print(f"|dE/dx| = (2 * {E0}) / (3 * {R0} * (1 - {x} / {R0})^(1/3))")
print(f"\nThe calculated energy loss per centimeter at {x} cm is: {dEdx:.3f} MeV/cm")
