import math

# Given parameters
# rho = 1000 kg/m^3 (density of water)
# g = 10 m/s^2 (acceleration due to gravity)
# H = 10 m (depth of the river)

g = 10
H = 10
rho = 1000

# The initial gauge pressure at the bottom of the river is P_static = rho * g * H
# When the water flows, the pressure is reduced. We want to find the speed v
# where the gauge pressure at the bottom becomes 0.
# The condition is: rho * g * H = (1/2) * rho * v^2
# This simplifies to: g * H = 0.5 * v^2
# Or: v = sqrt(2 * g * H)

print("To find the flow speed v, we solve the equation derived from Bernoulli's principle:")
print("Pressure reduction due to flow = Initial hydrostatic pressure")
print(f"(1/2) * rho * v^2 = rho * g * H")
print("We can simplify this to:")
print("v = sqrt(2 * g * H)")
print("\nSubstituting the given values:")
two = 2
v_squared = two * g * H
v = math.sqrt(v_squared)
print(f"v = sqrt({two} * {g} * {H})")
print(f"v = sqrt({v_squared})")
print(f"The flow speed v is {v} m/s.")
