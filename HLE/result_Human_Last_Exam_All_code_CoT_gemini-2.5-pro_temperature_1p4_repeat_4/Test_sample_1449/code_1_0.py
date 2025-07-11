import math

# --- Given constants ---
m = 0.1  # kg (mass of the block, 100g)
M = 10.0 # kg (mass of the wedge)
h = 2.0  # m (height of the wedge)
theta_deg = 30.0 # degrees (angle of the wedge)
mu = 0.5 # coefficient of friction
g = 10.0 # m/s^2 (acceleration due to gravity)

# Convert angle to radians for trigonometric functions
theta_rad = math.radians(theta_deg)
s = math.sin(theta_rad)
c = math.cos(theta_rad)
cos2theta = math.cos(2 * theta_rad)

print("--- Step 1: Calculate the relative acceleration (a_rel) ---")

# The formula for relative acceleration a_rel, derived from the equations of motion:
# a_rel = g * [M*(s - μ*c) - m*cos(2θ)*(s + μ*c)] / [M + m*s*(s + μ*c)]
numerator = g * (M * (s - mu * c) - m * cos2theta * (s + mu * c))
denominator = M + m * s * (s + mu * c)
a_rel = numerator / denominator

print(f"The acceleration of the block relative to the wedge is calculated as:")
print(f"a_rel = (g * (M*(sinθ - μ*cosθ) - m*cos(2θ)*(sinθ + μ*cosθ))) / (M + m*sinθ*(sinθ + μ*cosθ))")
print(f"a_rel = ({g} * ({M}*({s:.3f} - {mu}*{c:.3f}) - {m}*{cos2theta:.3f}*({s:.3f} + {mu}*{c:.3f}))) / ({M} + {m}*{s:.3f}*({s:.3f} + {mu}*{c:.3f}))")
print(f"a_rel = {numerator:.4f} / {denominator:.4f}")
print(f"a_rel = {a_rel:.4f} m/s²\n")


print("--- Step 2: Calculate the distance the block slides (L) ---")
# The distance is the length of the wedge's hypotenuse
L = h / s
print(f"The distance L is calculated using L = h / sin(θ)")
print(f"L = {h} / {s:.3f} = {L:.4f} m\n")

print("--- Step 3: Calculate the time taken (t) ---")
# Using the kinematic equation L = 0.5 * a_rel * t^2
# t = sqrt(2 * L / a_rel)
time = math.sqrt(2 * L / a_rel)

print(f"The time t is calculated using t = sqrt(2 * L / a_rel)")
print(f"t = sqrt(2 * {L:.4f} / {a_rel:.4f})")
print(f"t = {time:.4f} seconds")

# Final Answer
# print(f"\n<<< {time} >>>")