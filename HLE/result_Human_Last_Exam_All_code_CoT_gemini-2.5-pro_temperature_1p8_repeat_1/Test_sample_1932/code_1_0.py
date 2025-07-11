import math

# Description of the model and derivation
print("Analysis of the Hourglass Weight Change:")
print("------------------------------------------")
print("An hourglass's weight changes while running due to two main effects:")
print("1. Weight of airborne sand: This mass is in free-fall and not supported, making the hourglass lighter.")
print("2. Impact force of landing sand: This imparts a downward force, making the hourglass heavier.")
print("\nIn a simple model, these effects cancel. However, for granular sand, the landing is 'soft'.")
print("The impact force transmitted to the base is less than the weight of the airborne sand.")
print("Therefore, the hourglass becomes slightly LIGHTER.")
print("\nWe estimate this change by the weight of the airborne sand, using the maximum fall distance H for the 'largest possible effect'.")

# Parameters (converted to SI units)
d = 0.01  # m (1 cm)
h = 0.02  # m (2 cm)
H = 0.04  # m (4 cm)
rho = 1500  # kg/m^3
t = 60  # s (1 minute)
g = 9.8  # m/s^2

# The chosen formula from the options
print("\nThe formula for the change in weight is:")
print("ΔW = - (π * d² * h * ρ / (4 * t)) * √(2 * g * H)")
print("------------------------------------------")


# Break down the calculation
M_sand = rho * (math.pi * d**2 / 4) * h
m_dot = M_sand / t
v_impact = math.sqrt(2 * g * H)
delta_W = -m_dot * v_impact

# Print the equation with numerical values plugged in
print("\nPlugging in the numbers:")
# The format below shows each number in the equation as requested.
print(f"Mass flow rate (ṁ) = (π * {d}² * {h} * {rho}) / (4 * {t}) = {m_dot:.4e} kg/s")
print(f"Impact velocity (v) = √(2 * {g} * {H}) = {v_impact:.4f} m/s")
print(f"ΔW = - {m_dot:.4e} kg/s * {v_impact:.4f} m/s")

# Final result
print("\nFinal Result:")
print(f"The estimated change in weight is: {delta_W:.4e} Newtons.")
print(f"(A negative value means the hourglass is lighter while running).")
print(f"This is a weight change of approximately {abs(delta_W) * 1000:.4f} milliNewtons.")