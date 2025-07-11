import math

# Given parameters
m = 0.1  # kg
M = 10.0 # kg
theta_deg = 30.0 # degrees
h = 2.0  # meters
mu = 0.5
g = 10.0 # m/s^2

# Convert angle to radians for math functions
theta_rad = math.radians(theta_deg)

# --- Step 1: Calculate the sliding distance L ---
# The distance the block slides, L, is the length of the wedge's hypotenuse.
L = h / math.sin(theta_rad)
print("Step 1: Calculate the sliding distance L")
print(f"L = h / sin(θ)")
print(f"L = {h} / sin({theta_deg}°) = {L:.4f} m\n")

# --- Step 2: Calculate the relative acceleration a_rel ---
# The formula for the acceleration of the block relative to the wedge is derived from
# analyzing the forces on both the block and the wedge.
# a_rel = g * (sin(θ) - μ*cos(θ)) * (M+m) / (M + m*sin²(θ) - m*μ*sin(θ)cos(θ))
s_t = math.sin(theta_rad)
c_t = math.cos(theta_rad)

# Numerator of a_rel
numerator = g * (s_t - mu * c_t) * (M + m)
# Denominator of a_rel
denominator = M + m * s_t**2 - m * mu * s_t * c_t

a_rel = numerator / denominator

print("Step 2: Calculate the relative acceleration a_rel")
print("a_rel = g*(sin(θ) - μ*cos(θ))*(M+m) / (M + m*sin²(θ) - m*μ*sin(θ)cos(θ))")
print(f"a_rel = {g}*(sin({theta_deg}°) - {mu}*cos({theta_deg}°))*({M}+{m}) / ({M} + {m}*sin²({theta_deg}°) - {m}*{mu}*sin({theta_deg}°)*cos({theta_deg}°))")
print(f"a_rel = {g}*({s_t:.4f} - {mu}*{c_t:.4f})*({M+m}) / ({M} + {m}*{s_t**2:.4f} - {m}*{mu}*{s_t:.4f}*{c_t:.4f})")
print(f"a_rel = {numerator:.4f} / {denominator:.4f}")
print(f"a_rel = {a_rel:.4f} m/s²\n")

# --- Step 3: Calculate the time t to slide down ---
# Using the kinematic equation L = (1/2)*a_rel*t² for an object starting from rest.
# t = sqrt(2 * L / a_rel)
time_to_slide = math.sqrt(2 * L / a_rel)

print("Step 3: Calculate the time t")
print("t = sqrt(2 * L / a_rel)")
print(f"t = sqrt(2 * {L:.4f} / {a_rel:.4f})")
print(f"t = {time_to_slide:.4f} s")

# Final answer
final_answer = time_to_slide
print(f"\nFinal Answer: The exact time it takes for the block to slide down is {final_answer:.4f} seconds.")
print(f"<<<{final_answer}>>>")