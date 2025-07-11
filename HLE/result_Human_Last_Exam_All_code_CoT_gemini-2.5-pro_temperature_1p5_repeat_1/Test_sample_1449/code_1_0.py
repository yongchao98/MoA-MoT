import math

# --- Given Parameters ---
m = 100 / 1000  # Convert g to kg (0.1 kg)
M = 10.0        # kg
h = 2.0         # m
theta_deg = 30.0  # degrees
mu = 1.0 / 2.0    # coefficient of friction
g = 10.0        # m/s^2

# The initial velocity v0=1 m/s of the wedge (and block) does not affect
# the relative acceleration or the time to slide down, so it is not used in the calculation.

# --- Step-by-Step Calculation ---

# Convert angle to radians for trigonometric functions
theta_rad = math.radians(theta_deg)
s = math.sin(theta_rad)
c = math.cos(theta_rad)

# 1. Solve for the acceleration of the wedge (a_M) in the horizontal direction.
# This is derived by analyzing the forces between the block and the wedge.
# The horizontal force on the wedge is M*a_M = N*sin(theta) - f_k*cos(theta) where f_k = mu*N.
# However, N depends on a_M. The final combined equation for a_M is:
# a_M = (m*g*c*(mu*c - s)) / (M - m*s*(mu*c - s))
numerator_a_M = m * g * c * (mu * c - s)
denominator_a_M = M - m * s * (mu * c - s)
a_M = numerator_a_M / denominator_a_M

# 2. Solve for the acceleration of the block relative to the wedge (a_rel).
# In the non-inertial frame of the wedge, the acceleration down the incline is given by:
# a_rel = g*(sin(theta) - mu*cos(theta)) - a_M*(mu*sin(theta) + cos(theta))
term1_arel = g * (s - mu * c)
term2_arel_factor = mu * s + c
a_rel = term1_arel - a_M * term2_arel_factor

# 3. Calculate the distance (L) the block travels down the incline.
# L is the hypotenuse of the triangle with height h.
L = h / s

# 4. Use kinematics to find the time (t).
# The block starts from rest relative to the wedge.
# Equation: L = v0*t + 0.5*a*t^2 => L = 0.5 * a_rel * t^2
# Solve for t: t = sqrt(2 * L / a_rel)
time_squared = (2 * L) / a_rel
time = math.sqrt(time_squared)

# --- Output the results with equations ---
print("--- Final Calculation Breakdown ---")

print("\n1. First, calculate the distance the block slides (L):")
print(f"L = h / sin(theta) = {h:.1f} / sin({theta_deg}) = {h:.1f} / {s:.3f} = {L:.4f} m")

print("\n2. Next, calculate the acceleration of the wedge (a_M):")
print("a_M = (m*g*cos(th)*(mu*cos(th)-sin(th))) / (M - m*sin(th)*(mu*cos(th)-sin(th)))")
print(f"a_M = ({m:.1f}*{g:.1f}*{c:.3f}*({mu:.1f}*{c:.3f}-{s:.3f})) / ({M:.1f} - {m:.1f}*{s:.3f}*({mu:.1f}*{c:.3f}-{s:.3f}))")
print(f"a_M = {numerator_a_M:.4f} / {denominator_a_M:.4f} = {a_M:.4f} m/s^2")

print("\n3. Then, calculate the block's acceleration relative to the wedge (a_rel):")
print("a_rel = g*(sin(th) - mu*cos(th)) - a_M*(mu*sin(th) + cos(th))")
print(f"a_rel = {g:.1f}*({s:.3f} - {mu:.1f}*{c:.3f}) - ({a_M:.4f})*({mu:.1f}*{s:.3f} + {c:.3f})")
print(f"a_rel = {term1_arel:.4f} - ({a_M * term2_arel_factor:.4f}) = {a_rel:.4f} m/s^2")

print("\n4. Finally, calculate the total time (t) using kinematics:")
print("t = sqrt(2 * L / a_rel)")
print(f"t = sqrt(2 * {L:.4f} / {a_rel:.4f})")
print(f"t = sqrt({time_squared:.4f})")
print(f"\nt = {time:.4f} s")

# The final answer in the requested format
final_answer_value = f"{time:.4f}"
# The below print is just for final display and would not be part of the actual value.
# print(f"\nFinal answer formatted for submission: <<<{final_answer_value}>>>")