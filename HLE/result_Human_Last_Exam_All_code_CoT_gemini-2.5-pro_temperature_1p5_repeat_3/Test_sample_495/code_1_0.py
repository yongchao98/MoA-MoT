import math

# --- Initial values ---
E0 = 8.5  # Initial energy in MeV
R0 = 8.3  # Total range in cm
x = 4.0   # Distance from the source in cm

# --- Step 1: Calculate the constant 'a' from Geiger's rule R = a * E^(3/2) ---
# a = R0 / E0^(3/2)
a = R0 / (E0**1.5)
print(f"Step 1: Calculate the constant 'a'")
print(f"Using R = a * E^(3/2), with R0 = {R0} cm and E0 = {E0} MeV:")
print(f"a = {R0} / {E0}**1.5 = {a:.4f} cm/MeV^(3/2)\n")

# --- Step 2: Find the remaining range and the energy at distance x ---
# The remaining range Rx at distance x is R0 - x
Rx = R0 - x
# The energy Ex at this point is found by rearranging Geiger's rule: Ex = (Rx / a)^(2/3)
Ex = (Rx / a)**(2.0/3.0)
print(f"Step 2: Find the energy at x = {x} cm")
print(f"Remaining range, Rx = {R0} cm - {x} cm = {Rx} cm")
print(f"Energy at {x} cm, Ex = ({Rx} / {a:.4f})^(2/3) = {Ex:.4f} MeV\n")

# --- Step 3: Calculate the energy loss per centimeter, |dE/dx| ---
# The energy loss is given by |dE/dx| = 1 / (a * (3/2) * sqrt(E))
energy_loss_per_cm = 1 / (a * 1.5 * math.sqrt(Ex))
print(f"Step 3: Calculate the energy loss per cm")
print(f"Using the formula |dE/dx| = 1 / (a * 1.5 * E_x^0.5)")
print(f"Energy Loss = 1 / ({a:.4f} * 1.5 * {Ex:.4f}^0.5)")
print(f"Energy Loss = 1 / ({a:.4f} * 1.5 * {math.sqrt(Ex):.4f})")
print(f"Energy Loss = 1 / ({a * 1.5 * math.sqrt(Ex):.4f})")
print(f"\nThe calculated energy loss per centimetre at {x} cm is: {energy_loss_per_cm:.3f} MeV/cm")

# --- Final Answer ---
final_answer = energy_loss_per_cm
# The output format is just the number, so we will not print the final string here
# but the calculation above leads to the final number.