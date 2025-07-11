import math

# --- Given Parameters ---
# The problem states the river is 10 meters deep.
H = 10  # meters

# The problem uses g = 10 m/s^2 in its example calculation.
g = 10  # m/s^2

# The density of water is typically 1000 kg/m^3.
rho = 1000 # kg/m^3

# --- Explanation and Calculation ---

# According to Bernoulli's principle, the pressure in a moving fluid decreases as its speed increases.
# The initial gauge pressure at the bottom of the stationary river is the hydrostatic pressure P.
# P = ρ * g * H
initial_pressure = rho * g * H

print(f"The initial pressure at the bottom of the river due to the water column is:")
print(f"P = ρ * g * H = {rho} * {g} * {H} = {int(initial_pressure)} N/m^2\n")

# When the water flows at speed v, the pressure drops by an amount equal to the dynamic pressure, (1/2) * ρ * v^2.
# We want to find the speed 'v' at which the pressure at the bottom becomes zero.
# This occurs when the dynamic pressure drop is equal to the initial hydrostatic pressure.
print("To find the speed 'v' that makes the bottom pressure zero, we set the initial pressure equal to the dynamic pressure drop:")
print("ρ * g * H = (1/2) * ρ * v^2\n")

# The density ρ cancels from both sides of the equation.
print("Simplifying the equation, we get:")
print("g * H = (1/2) * v^2")
print("v^2 = 2 * g * H")
print("v = sqrt(2 * g * H)\n")

# Now, we substitute the given values into the equation to find v.
print("Plugging in the values for g and H:")
print(f"v = sqrt(2 * {g} * {H})")

# Calculate the final velocity
v = math.sqrt(2 * g * H)

print(f"v = sqrt({2 * g * H})")
print(f"v ≈ {v:.2f} m/s")

# --- Final Answer ---
# The problem asks for the final numerical answer in a specific format.
print(f"\n<<<{v:.2f}>>>")