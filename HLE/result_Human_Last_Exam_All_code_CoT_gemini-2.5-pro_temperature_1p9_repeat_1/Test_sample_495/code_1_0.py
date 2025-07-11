import math

# --- Given values ---
# Initial energy in MeV
E0 = 8.5
# Total range in cm
R0 = 8.3
# Distance from the source in cm
x = 4.0

# --- Explanation of the Formula ---
print("The calculation for the energy loss per centimeter, |dE/dx|, is based on the empirical range-energy relationship R = a * E^(3/2).")
print("From this, the rate of energy loss at a distance x from the source can be expressed as:")
print("|dE/dx| = (2 * E0) / (3 * R0^(2/3) * (R0 - x)^(1/3))\n")

# --- Displaying the Equation with Substituted Values ---
print("Substituting the given values into the formula:")
# Note: Here we explicitly show each number that will be used in the calculation
print(f"|dE/dx| = (2 * {E0}) / (3 * {R0}^(2/3) * ({R0} - {x})^(1/3))")
# Calculate the remaining range
remaining_range = R0 - x
print(f"|dE/dx| = ({2 * E0}) / (3 * {R0}^(2/3) * {remaining_range}^(1/3))\n")


# --- Calculation ---
# Calculate the energy loss per centimeter (|dE/dx|)
denominator = 3 * (R0**(2/3)) * (remaining_range**(1/3))
dEdx = (2 * E0) / denominator


# --- Final Result ---
print(f"The energy loss per centimetre for the α-particles at a distance of {x} cm from the source is: {dEdx:.4f} MeV/cm.")
