import math

# Plan:
# 1. Define all given variables from the problem statement and the table.
# 2. State the main assumption for solving the problem: applying Bernoulli's equation
#    between the tank surface and the manometer.
# 3. Calculate intermediate values like pipe diameter, area, and fluid velocity.
# 4. Calculate all components of the Bernoulli equation: pressure head, kinetic head,
#    elevation head, and head loss (friction + entrance).
# 5. Sum these components to find the total pump head (h_pump).
# 6. Calculate the work of the pump (power in Watts) using the pump head.
# 7. Print each step of the calculation, showing the formulas and the numbers used.

# Step 1: Define variables from the problem
rho = 997  # density of water in kg/m^3
P1 = 1.05e5  # ambient pressure (at tank surface) in N/m^2
P_p = 1.33e5  # manometer absolute pressure in N/m^2
Q = 2.86e-3  # volume flow rate in m^3/s
r = 15.5 / 1000  # pipe radius in meters
L = 14.9  # pipe length in meters
f = 0.004  # Darcy friction factor
K_c = 0.4  # entrance loss coefficient (shrinkage)
g = 9.81  # acceleration due to gravity in m/s^2

# Step 2: State assumptions based on problem interpretation
# Assumption: The elevation difference between the water surface in the tank (z1)
# and the pump/manometer (zp) is 2 m, based on "The amount of water in the tank is 2 meters".
z1_minus_zp = 2.0  # m

# Step 3: Calculate pipe geometric properties and fluid velocity
D = 2 * r  # pipe diameter in m
A = math.pi * r**2  # pipe cross-sectional area in m^2
v = Q / A  # velocity in m/s

# Step 4: Calculate the terms for the extended Bernoulli equation
# Equation: h_pump = (P_p - P1)/(ρg) + v^2/(2g) + (z_p - z1) + h_loss
pressure_head = (P_p - P1) / (rho * g)
kinetic_head = v**2 / (2 * g)
elevation_head_diff = -z1_minus_zp
friction_loss = f * (L / D) * kinetic_head
entrance_loss = K_c * kinetic_head
total_head_loss = friction_loss + entrance_loss

# Step 5: Calculate the pump head (h_pump)
h_pump = pressure_head + kinetic_head + elevation_head_diff + total_head_loss

# Step 6: Calculate the work of the pump (Power in Watts)
work_pump = h_pump * rho * g * Q

# Step 7: Print the results step-by-step
print("--- Calculation Steps for Pump Work ---")
print(f"1. Pipe Diameter D = 2 * {r:.4f} m = {D:.4f} m")
print(f"2. Pipe Area A = \u03C0 * ({r:.4f} m)\u00b2 = {A:.6f} m\u00b2")
print(f"3. Water Velocity v = {Q:.5f} m\u00b3/s / {A:.6f} m\u00b2 = {v:.3f} m/s")
print("\n--- Applying Bernoulli's Equation to find Pump Head (h_pump) ---")
print("Equation: h_pump = (P_p - P1)/(\u03C1*g) + v\u00b2/(2*g) + (z_p - z1) + h_loss")

print("\n- Calculating Components:")
print(f"  - Pressure Head = ({P_p:.0f} - {P1:.0f}) / ({rho} * {g}) = {pressure_head:.3f} m")
print(f"  - Kinetic Head = {v:.3f}\u00b2 / (2 * {g}) = {kinetic_head:.3f} m")
print(f"  - Elevation Head (z_p - z1) = {elevation_head_diff:.1f} m")
print(f"  - Head Loss (h_loss = friction + entrance):")
print(f"    - Friction Loss = {f} * ({L}/{D:.4f}) * {kinetic_head:.3f} m = {friction_loss:.3f} m")
print(f"    - Entrance Loss = {K_c} * {kinetic_head:.3f} m = {entrance_loss:.3f} m")
print(f"    - Total Head Loss h_loss = {friction_loss:.3f} m + {entrance_loss:.3f} m = {total_head_loss:.3f} m")

print("\n- Calculating Final Pump Head:")
print(f"h_pump = {pressure_head:.3f} m + {kinetic_head:.3f} m + ({elevation_head_diff:.1f} m) + {total_head_loss:.3f} m")
print(f"h_pump = {h_pump:.3f} m")

print("\n--- Final Pump Work (Power) Calculation ---")
print("Equation: W_pump = h_pump * \u03C1 * g * Q")
print(f"W_pump = {h_pump:.3f} m * {rho} kg/m\u00b3 * {g} m/s\u00b2 * {Q:.5f} m\u00b3/s")
print(f"W_pump = {work_pump:.2f} Watts")
<<<92.15>>>