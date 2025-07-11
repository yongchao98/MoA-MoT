import math

# --- Given Information ---
# Initial energy of the alpha-particle in MeV
E0 = 8.5
# Total range of the alpha-particle in cm
R0 = 8.3
# Distance from the source where we want to find the energy loss in cm
x = 4.0

# --- Step 1: Calculate the constant 'k' ---
# The range-energy relationship is R = k * E^(3/2)
# We can find k using the initial conditions: k = R0 / E0^(3/2)
k = R0 / (E0 ** 1.5)

# --- Step 2: Find the energy of the particle at distance x ---
# First, calculate the remaining range of the particle at distance x
remaining_range = R0 - x
# Now, use the range-energy formula to find the energy (Ex) at this point
# remaining_range = k * Ex^(3/2)  =>  Ex = (remaining_range / k)^(2/3)
Ex = (remaining_range / k) ** (2/3)

# --- Step 3: Calculate the energy loss per centimeter ---
# The formula for energy loss per cm (-dE/dx) derived from the range-energy relation is:
# -dE/dx = 1 / (1.5 * k * Ex^(1/2))
energy_loss_per_cm = 1 / (1.5 * k * (Ex ** 0.5))

# --- Output the results step-by-step ---
print("Calculating the energy loss per centimetre for an alpha-particle.")
print(f"Initial Energy (E0): {E0} MeV")
print(f"Total Range (R0): {R0} cm")
print(f"Distance from source (x): {x} cm")
print("-" * 30)

print(f"Step 1: Calculate the proportionality constant 'k'.")
print(f"k = R0 / E0^(3/2) = {R0} / {E0}^(1.5) = {k:.4f} cm/MeV^(3/2)")
print("-" * 30)

print(f"Step 2: Find the energy at x = {x} cm.")
print(f"The remaining range is {R0} - {x} = {remaining_range:.1f} cm.")
print("The energy (Ex) at this point is calculated as (remaining_range / k)^(2/3).")
print(f"Ex = ({remaining_range:.1f} / {k:.4f})^(2/3) = {Ex:.4f} MeV")
print("-" * 30)

print(f"Step 3: Calculate the energy loss per cm (-dE/dx).")
print("The formula is: -dE/dx = 1 / (1.5 * k * Ex^(1/2))")
print(f"-dE/dx = 1 / (1.5 * {k:.4f} * {Ex:.4f}^(0.5))")
final_equation_denominator = 1.5 * k * (Ex ** 0.5)
print(f"-dE/dx = 1 / {final_equation_denominator:.4f}")
print("-" * 30)
print(f"Final Answer: The energy loss per centimetre at {x} cm is {energy_loss_per_cm:.2f} MeV/cm.")
