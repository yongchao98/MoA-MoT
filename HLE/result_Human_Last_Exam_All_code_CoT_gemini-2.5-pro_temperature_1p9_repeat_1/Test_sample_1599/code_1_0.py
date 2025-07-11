# Initial values
u = 1.5  # m/s
t0 = 0.0
t1 = 4.0
a1 = -0.15  # m/s^2
t3 = 15.0
a2 = 0.25  # m/s^2
t5 = 23.0
t6 = 40.0

# --- Step 1: Calculate the man's position and velocity up to t5=23s ---

# Position and velocity at t1
p_m_t1 = u * (t1 - t0)
v_m_t1 = u
print(f"Man's state at t1={t1}s: position = {p_m_t1} m, velocity = {v_m_t1} m/s")

# Calculate t2, when the man comes to a stop
# v_m(t2) = v_m(t1) + a1 * (t2 - t1) = 0
t2 = t1 - v_m_t1 / a1
dt_12 = t2 - t1
print(f"Man comes to a stop at t2={t2}s")

# Position at t2
p_m_t2 = p_m_t1 + v_m_t1 * dt_12 + 0.5 * a1 * dt_12**2
v_m_t2 = 0
print(f"Man's position at t2={t2}s: position = {p_m_t2} m")

# Position and velocity at t3
p_m_t3 = p_m_t2
v_m_t3 = v_m_t2
print(f"Man's state at t3={t3}s: position = {p_m_t3} m, velocity = {v_m_t3} m/s")

# Calculate t4, when the man gets back to speed u
# v_m(t4) = v_m(t3) + a2 * (t4 - t3) = u
t4 = t3 + (u - v_m_t3) / a2
dt_34 = t4 - t3
print(f"Man reaches speed u={u}m/s at t4={t4}s")

# Position at t4
p_m_t4 = p_m_t3 + v_m_t3 * dt_34 + 0.5 * a2 * dt_34**2
v_m_t4 = u
print(f"Man's state at t4={t4}s: position = {p_m_t4} m, velocity = {v_m_t4} m/s")

# Position and velocity at t5
dt_45 = t5 - t4
p_m_t5 = p_m_t4 + v_m_t4 * dt_45
v_m_t5 = v_m_t4
print(f"Man's state at t5={t5}s: position = {p_m_t5} m, velocity = {v_m_t5} m/s")


# --- Step 2: Set up the equation for the man's final position at t6 ---

dt_56 = t6 - t5
# Final position y_meet = p_m_t5 + v_m_t5 * dt_56 + 0.5 * a3 * dt_56**2
# This simplifies to y_meet = (p_m_t5 + v_m_t5 * dt_56) + (0.5 * dt_56**2) * a3
const_term_man = p_m_t5 + v_m_t5 * dt_56
coeff_a3 = 0.5 * dt_56**2

print(f"\nFinal meeting position for the man is y_meet = {const_term_man:.1f} + {coeff_a3:.1f} * a3")


# --- Step 3: Use bird's motion to find another expression for y_meet ---

# Based on the reasoning, the planned rendezvous time T_planned was also t=40s.
T_planned = 40.0
# The planned meeting point's y-coordinate is based on the man's planned (constant velocity) motion from t4.
y_meet_planned = p_m_t4 + v_m_t4 * (T_planned - t4)
print(f"Planned meeting y-coordinate was {y_meet_planned} m")

# The bird's y-position at the time of the gust is the same as the planned meeting y-coordinate.
y_b_t5 = y_meet_planned

# The northward speed of the bird, vw, is assumed to be the man's characteristic speed u.
v_w = u
print(f"Assuming wind gust gives bird northward velocity v_w = {v_w} m/s")

# The bird's final y-position at t6 is:
y_meet_bird = y_b_t5 + v_w * dt_56
print(f"Final meeting position from bird's motion is y_meet = {y_meet_bird:.1f} m")

# --- Step 4: Solve for the unknown acceleration a3 ---
# Equating the two expressions for y_meet:
# const_term_man + coeff_a3 * a3 = y_meet_bird
# Note: const_term_man is 21 + 1.5*17 = 46.5, and y_meet_planned is 18 + 1.5*19 = 46.5. They are equal.
# So const_term_man + coeff_a3 * a3 = y_b_t5 + v_w * dt_56 simplifies to:
# coeff_a3 * a3 = v_w * dt_56
a3 = (v_w * dt_56) / coeff_a3
print(f"\nFinal equation: {coeff_a3:.1f} * a3 = {v_w} * {dt_56}")
print(f"Solving for a3: a3 = ({v_w} * {dt_56}) / {coeff_a3:.1f}")

# Numerator and denominator for fraction simplification
numerator = 17 * 1.5 * 2 # Multiply by 2 to deal with 0.5 in denominator
denominator = 17**2
common_factor = 17
num_simple = int(numerator / common_factor)
den_simple = int(denominator / common_factor)


print(f"The man's constant acceleration during this final leg is {a3} m/s^2.")
print(f"As a fraction, a3 = {num_simple}/{den_simple} m/s^2.")

# <<< [final_answer] >>>
# Return the numerical value of a3
final_answer = a3
print(f'<<<_A>>> {final_answer}')