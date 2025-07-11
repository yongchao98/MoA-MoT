import math

# Description of the model and derivation of the formula.
# The change in weight (Delta_W) of a running hourglass can be estimated by considering the dominant physical effects.
# While an idealized model suggests the weight remains unchanged, a more realistic analysis points to the
# impact force of the falling sand stream as the largest effect. This makes the hourglass effectively heavier.
#
# The change in weight is thus approximated by the impact force:
# Delta_W ≈ F_impact = mass_flow_rate * impact_velocity
#
# 1. Mass flow rate (m_dot): The total mass of sand (M = rho * V) falls in time t.
#    The volume V of the sand is that of a cylinder with diameter d and height h: V = (pi * d^2 / 4) * h.
#    So, m_dot = M / t = (pi * d^2 * h * rho) / (4 * t).
#
# 2. Impact velocity (v_f): The sand falls from the top chamber to the bottom one. The fall distance
#    changes, but we can use the characteristic maximum fall distance, which is the chamber height H.
#    Assuming free-fall, the velocity upon impact is v_f = sqrt(2 * g * H).
#
# Combining these gives the final expression for the estimated weight change:
# Delta_W = (pi * d^2 * h * rho / (4 * t)) * sqrt(2 * g * H)

# Given parameters with units converted to SI
d = 0.01      # diameter in meters
h = 0.02      # sand pile height in meters
H = 0.04      # chamber height in meters
rho = 1500    # density in kg/m^3
t = 60        # time in seconds
g = 9.8       # acceleration due to gravity in m/s^2

# Calculation steps
# Term 1: Mass flow rate
mass_flow_rate_numerator = math.pi * d**2 * h * rho
mass_flow_rate_denominator = 4 * t
mass_flow_rate = mass_flow_rate_numerator / mass_flow_rate_denominator

# Term 2: Impact velocity
impact_velocity = math.sqrt(2 * g * H)

# Final calculation for Delta_W
delta_W = mass_flow_rate * impact_velocity

# Output the step-by-step calculation with numbers
print("Final Expression: Delta_W = (pi * d^2 * h * rho / (4 * t)) * sqrt(2 * g * H)\n")

print(f"Step 1: Calculate the mass flow rate (m_dot)")
print(f"m_dot = (pi * {d}^2 * {h} * {rho}) / (4 * {t})")
print(f"m_dot = ({mass_flow_rate_numerator:.4e}) / ({mass_flow_rate_denominator})")
print(f"m_dot = {mass_flow_rate:.4e} kg/s\n")

print(f"Step 2: Calculate the impact velocity (v_f)")
print(f"v_f = sqrt(2 * {g} * {H})")
print(f"v_f = sqrt({2 * g * H:.4f})")
print(f"v_f = {impact_velocity:.4f} m/s\n")

print(f"Step 3: Calculate the change in weight (Delta_W)")
print(f"Delta_W = m_dot * v_f")
print(f"Delta_W = {mass_flow_rate:.4e} kg/s * {impact_velocity:.4f} m/s")
print(f"Delta_W = {delta_W:.4e} N")