import math

# --- Given Data ---
# Side of the equilateral triangle (m)
a = 1.2 * 10**10
# Tangential velocity (m/s), converted from km/s
v = 125 * 1000
# Gravitational constant (kg^-1 m^3 s^-2)
G = 6.67 * 10**-11
# Solar mass (kg)
M_solar = 1.99 * 10**30

# --- Calculation ---
# As derived in the plan, the formula for the mass 'm' of a single component is:
# m = a * v^2 / G

# Calculate the mass in kilograms
m_kg = (a * v**2) / G

# Convert the mass to solar masses
m_solar = m_kg / M_solar

# --- Output ---
print("The formula to find the mass of a single star component is m = a * v^2 / G.")
print("Plugging in the given values, the equation is:")
print(f"m = ({a}) * ({v})^2 / ({G})")

print(f"\nCalculated mass in kilograms: {m_kg:.3e} kg")

print("\nConverting to solar masses by dividing by the solar mass:")
print(f"Mass in solar masses = {m_kg:.3e} kg / {M_solar:.3e} kg")
print(f"Mass in solar masses = {m_solar:.3f}")

# Final answer rounded to one decimal place
final_answer = round(m_solar, 1)

print(f"\nRounded to one decimal place, the mass of a single component is {final_answer} solar masses.")