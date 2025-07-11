import math

# --- Given Data and Constants ---
# Side length of the equilateral triangle in meters
s = 1.2 * 10**10
# Tangential velocity in km/s
v_km_s = 125
# Mass of the sun in kg
M_sun = 1.99 * 10**30
# Gravitational constant in m^3 kg^-1 s^-2
G = 6.67 * 10**-11

# --- Calculations ---
# Convert velocity from km/s to m/s for use in SI units
v_m_s = v_km_s * 1000

# As derived in the plan, the mass 'm' of a single component can be found
# by equating the net gravitational force on a star to the centripetal force
# required for its orbit. This yields the formula: m = (v^2 * s) / G.
# We use the given instantaneous velocity 'v' and side length 's'.

# Calculate the mass of a single star component in kilograms
m_kg = (v_m_s**2 * s) / G

# Convert the mass from kilograms to solar masses
m_solar = m_kg / M_sun

# --- Output the results step-by-step ---
print("Step 1: The equation to find the mass 'm' of a single star is derived from F_gravity = F_centripetal.")
print("The final equation is: m = (v^2 * s) / G\n")

print("Step 2: We substitute the given values into the equation.")
# We output each number used in the final equation.
print(f"m = ({v_m_s:.0f} m/s)^2 * ({s:.1e} m) / ({G:.2e} m^3 kg^-1 s^-2)\n")

print("Step 3: The calculated mass of a single component in kilograms is:")
print(f"m = {m_kg:.4e} kg\n")

print("Step 4: This mass is converted to solar masses by dividing by the Sun's mass.")
print(f"m (in solar masses) = {m_kg:.4e} kg / {M_sun:.2e} kg")
print(f"m = {m_solar:.4f} solar masses\n")

print("Step 5: The final answer, accurate to one decimal place, is:")
# Round the final result to one decimal place as requested
final_answer = round(m_solar, 1)
print(f"The mass of a single component is {final_answer} solar masses.")
