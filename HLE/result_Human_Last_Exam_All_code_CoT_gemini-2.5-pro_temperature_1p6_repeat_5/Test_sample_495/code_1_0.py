import math

# Initial given values
E0 = 8.5  # Initial energy in MeV
R0 = 8.3  # Total range in cm
x = 4.0   # Distance from the source in cm

print("Problem: Calculate the energy loss per centimeter for an alpha-particle at a distance of 4 cm from its source.")
print(f"Initial Energy (E0): {E0} MeV")
print(f"Total Range (R0): {R0} cm")
print(f"Distance from source (x): {x} cm")
print("-" * 50)

# Step 1: Calculate the constant k from Geiger's rule: R = k * E^(3/2)
# k = R0 / E0^(3/2)
k = R0 / (E0**1.5)
print("Step 1: Calculate the constant 'k' for Geiger's rule.")
print(f"Equation: k = R0 / E0^(3/2)")
print(f"Calculation: k = {R0:.1f} / {E0:.1f}^1.5 = {k:.4f} cm/MeV^1.5\n")

# Step 2: Calculate the remaining range (R) at distance x
R_at_x = R0 - x
print("Step 2: Calculate the remaining range 'R' at 4 cm.")
print(f"Equation: R = R0 - x")
print(f"Calculation: R = {R0:.1f} - {x:.1f} = {R_at_x:.1f} cm\n")

# Step 3: Calculate the energy (E) at that point using the remaining range R
# E = (R / k)^(2/3)
E_at_x = (R_at_x / k)**(2.0/3.0)
print("Step 3: Calculate the energy 'E' at 4 cm.")
print(f"Equation: E = (R / k)^(2/3)")
print(f"Calculation: E = ({R_at_x:.1f} / {k:.4f})^(2/3) = {E_at_x:.4f} MeV\n")

# Step 4: Calculate the energy loss per centimeter, |dE/dx|
# The derivative of the range-energy relation gives: |dE/dx| = (2/3) * (E / R)
dEdx = (2.0/3.0) * (E_at_x / R_at_x)
print("Step 4: Calculate the energy loss per centimeter, |dE/dx|.")
print(f"Equation: |dE/dx| = (2/3) * (E / R)")
print(f"Calculation: |dE/dx| = (2/3) * ({E_at_x:.4f} / {R_at_x:.1f}) = {dEdx:.4f} MeV/cm\n")

print("-" * 50)
print(f"Final Answer: The energy loss per centimetre at a distance of {x} cm is {dEdx:.3f} MeV/cm.")

# Format the final answer as requested
final_answer = round(dEdx, 3)
print(f"<<<{final_answer}>>>")