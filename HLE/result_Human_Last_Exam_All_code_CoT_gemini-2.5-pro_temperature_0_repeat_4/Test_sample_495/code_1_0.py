import math

# Initial conditions
E0 = 8.5  # Initial energy in MeV
R0 = 8.3  # Total range in cm
x = 4.0   # Distance from the source in cm

# Step 1: Calculate the constant k for standard air using R = k * E^(3/2)
# k = R0 / E0^(3/2)
k = R0 / (E0**1.5)

# Step 2: Find the energy of the particle at distance x = 4 cm.
# First, calculate the remaining range at distance x.
Rx = R0 - x
# Now, find the energy Ex corresponding to this remaining range using E = (R/k)^(2/3).
Ex = (Rx / k)**(2.0/3.0)

# Step 3: Calculate the energy loss per centimeter, |dE/dx|, at distance x.
# The formula is |dE/dx| = 1 / (1.5 * k * Ex^0.5)
dEdx = 1 / (1.5 * k * (Ex**0.5))

# Step 4: Print the final calculation and result.
print(f"The energy loss per centimetre for these α-particles at a distance of {x} cm is calculated as follows:")
print(f"Energy loss |dE/dx| = 1 / (1.5 * k * E(x)^0.5)")
print(f"|dE/dx| = 1 / (1.5 * {k:.4f} * {Ex:.4f}^0.5)")
print(f"|dE/dx| = 1 / (1.5 * {k:.4f} * {math.sqrt(Ex):.4f})")
print(f"|dE/dx| = {dEdx:.4f} MeV/cm")

# Final answer in the required format
final_answer = round(dEdx, 2)
# The final answer is printed below for clarity, but the required format is at the end.
# print(f"\nThe calculated energy loss is approximately {final_answer} MeV/cm.")