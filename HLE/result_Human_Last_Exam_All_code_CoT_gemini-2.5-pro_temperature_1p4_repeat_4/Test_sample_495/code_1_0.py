import math

# Step 1: Define initial values from the problem statement.
E0 = 8.5  # Initial energy in MeV
R0 = 8.3  # Total range in cm
x = 4.0   # Distance from the source in cm

# Step 2: Calculate the remaining range of the alpha particle.
# The remaining range (R_prime) is the total range minus the distance traveled.
R_prime = R0 - x

# Step 3: Calculate the energy loss per centimetre.
# The formula is derived from the empirical range-energy relationship R = cE^(3/2).
# Energy loss per cm (-dE/dx) = (2 * E0) / (3 * R0^(2/3) * R_prime^(1/3))
energy_loss = (2 * E0) / (3 * (R0**(2/3)) * (R_prime**(1/3)))

# Step 4: Print the explanation and the calculation.
print("To find the energy loss per centimetre, we use the formula derived from the range-energy relationship:")
print("-dE/dx = (2 * Initial_Energy) / (3 * Total_Range^(2/3) * Remaining_Range^(1/3))\n")

print("Given values:")
print(f"Initial Energy (E0): {E0} MeV")
print(f"Total Range (R0): {R0} cm")
print(f"Distance from source (x): {x} cm")
print(f"Remaining Range = {R0} cm - {x} cm = {R_prime} cm\n")

print("Plugging these values into the formula:")
print(f"-dE/dx = (2 * {E0}) / (3 * {R0}^(2/3) * {R_prime}^(1/3))")

# Step 5: Print the final result.
print(f"\nThe calculated energy loss per centimetre at {x} cm is: {energy_loss:.3f} MeV/cm")