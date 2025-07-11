import math

# Given parameters
d = 0.01  # m
h = 0.02  # m
H = 0.04  # m
rho = 1500  # kg/m^3
t = 60    # s
g = 9.8   # m/s^2

# The chosen expression for the weight change is from answer choice A.
# It is formulated as the product of mass flow rate (m_dot) and an impact velocity (v).
# delta_W = m_dot * v

# Calculate the mass flow rate (m_dot)
# m_dot = (total mass) / (total time)
# total mass = density * volume = rho * (pi * (d/2)^2 * h) = rho * pi * d^2 * h / 4
m_dot = (math.pi * d**2 * h * rho) / (4 * t)

# Calculate the characteristic impact velocity (v)
# v = sqrt(2 * g * L_eff) where L_eff is the effective fall height.
# For answer choice A, the effective fall height is taken as (H - h).
v = math.sqrt(2 * g * (H - h))

# Calculate the change in weight (delta_W)
delta_W = m_dot * v

# The problem asks for the expression and its components to be clear.
# Let's print the components of the final formula.
# The formula is: delta_W = (pi * d^2 * h * rho) / (4 * t) * sqrt(2 * g * (H - h))

print("The formula for the weight change (Delta W) corresponds to answer choice A:")
print(f"Delta W = (pi * d^2 * h * rho / (4 * t)) * sqrt(2 * g * (H - h))")
print("\nLet's calculate the value based on the given parameters:")
print("---------------------------------------------------------")
print(f"d = {d} m")
print(f"h = {h} m")
print(f"H = {H} m")
print(f"rho = {rho} kg/m^3")
print(f"t = {t} s")
print(f"g = {g} m/s^2")
print("---------------------------------------------------------")
print("The final equation with numerical values is:")
print(f"Delta W = (3.14159 * {d}^2 * {h} * {rho} / (4 * {t})) * sqrt(2 * {g} * ({H} - {h}))")
print("\nResult:")
print(f"The estimated change in weight is: {delta_W:.6e} Newtons")
