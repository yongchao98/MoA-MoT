import math

# --- Problem Data ---
# Side of the equilateral triangle (meters)
L = 1.2 * 10**10
# Tangential velocity (m/s), converted from km/s
v = 125 * 1000
# Gravitational constant (m^3 kg^-1 s^-2)
G = 6.67 * 10**-11
# Solar mass (kg)
M_solar = 1.99 * 10**30

# --- Calculation ---
# The net gravitational force on one star is F_net = sqrt(3) * G * m^2 / L^2.
# The centripetal force is F_c = m * v^2 / r, where r = L / sqrt(3) is the orbit radius.
# Equating F_net = F_c and solving for m yields the formula: m = (v^2 * L) / G

# Calculate the mass of a single component in kilograms
mass_kg = (v**2 * L) / G

# Convert the mass from kilograms to solar masses
mass_in_solar_masses = mass_kg / M_solar

# --- Output ---
print("To find the mass 'm' of a single star, we equate the net gravitational force with the centripetal force.")
print("The derived formula for the mass is: m = (v^2 * L) / G")
print("\nSubstituting the given values into the equation:")
print(f"m (kg) = ({v:.2e}^2 * {L:.2e}) / {G:.2e}")
mass_calculation_str = f"m (kg) = (({v**2:.2e}) * {L:.2e}) / {G:.2e}"
print(mass_calculation_str)
result_str = f"m (kg) = {(v**2 * L):.2e} / {G:.2e}"
print(result_str)
print(f"m = {mass_kg:.4e} kg")

print("\nConverting the mass to solar masses:")
print(f"Mass (solar masses) = Mass (kg) / Solar Mass (kg)")
print(f"Mass (solar masses) = {mass_kg:.4e} / {M_solar:.2e}")
print(f"\nThe mass of a single component is {mass_in_solar_masses:.1f} solar masses.")

# Final answer in the required format
final_answer = round(mass_in_solar_masses, 1)
print(f"\n<<<{final_answer}>>>")