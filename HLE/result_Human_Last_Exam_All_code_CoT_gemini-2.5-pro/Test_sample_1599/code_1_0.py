import math

# Step 1 & 2: Man's motion
# Given values
u = 1.5  # m/s
t0 = 0
t1 = 4
a1 = -0.15 # m/s^2
t3 = 15
a2 = 0.25 # m/s^2
t5 = 23
t6 = 40

# Man's position and velocity calculations
# At t0=0, y_m(0) = 0, u_m(0) = u
# Segment 1: t0 to t1 (constant velocity)
y_m_t1 = 0 + u * (t1 - t0)
u_m_t1 = u
print(f"Man's position at t1=4s: y_m({t1}) = {0} + {u} * ({t1} - {t0}) = {y_m_t1} m")

# Segment 2: t1 to t2 (deceleration to stop)
# v_f = v_i + a*t => 0 = u + a1*(t2-t1) => t2 = t1 - u/a1
t2 = t1 - u / a1
y_m_t2 = y_m_t1 + u_m_t1 * (t2 - t1) + 0.5 * a1 * (t2 - t1)**2
u_m_t2 = 0
print(f"Time man stops: t2 = {t1} - {u} / {a1} = {t2} s")
print(f"Man's position at t2={t2}s: y_m({t2}) = {y_m_t1} + {u_m_t1} * ({t2} - {t1}) + 0.5 * {a1} * ({t2} - {t1})^2 = {y_m_t2} m")

# Segment 3: t2 to t3 (at rest)
y_m_t3 = y_m_t2
u_m_t3 = u_m_t2
print(f"Man's position at t3={t3}s: y_m({t3}) = {y_m_t3} m (at rest)")

# Segment 4: t3 to t4 (acceleration to speed u)
# v_f = v_i + a*t => u = 0 + a2*(t4-t3) => t4 = t3 + u/a2
t4 = t3 + u / a2
y_m_t4 = y_m_t3 + u_m_t3 * (t4 - t3) + 0.5 * a2 * (t4 - t3)**2
u_m_t4 = u
print(f"Time man reaches speed u again: t4 = {t3} + {u} / {a2} = {t4} s")
print(f"Man's position at t4={t4}s: y_m({t4}) = {y_m_t3} + {u_m_t3} * ({t4} - {t3}) + 0.5 * {a2} * ({t4} - {t3})^2 = {y_m_t4} m")

# Segment 5: t4 to t5 (constant velocity)
y_m_t5 = y_m_t4 + u_m_t4 * (t5 - t4)
u_m_t5 = u_m_t4
print(f"Man's position at t5={t5}s: y_m({t5}) = {y_m_t4} + {u_m_t4} * ({t5} - {t4}) = {y_m_t5} m")

# Step 3, 4, 5: Connecting planned and actual motion
# The problem's structure implies that the planned rendezvous was set for the same time as the actual rendezvous, t6 = 40s.
# Duration of bird's planned final leg: T = t6 - t4
T = t6 - t4
print(f"Duration of bird's planned final leg: T = {t6} - {t4} = {T} s")

# The planned rendezvous y-coordinate is determined by the man's constant velocity travel from t4.
y_planned_rendezvous = y_m_t4 + u_m_t4 * T
print(f"Planned rendezvous y-coordinate: y_planned = {y_m_t4} + {u_m_t4} * {T} = {y_planned_rendezvous} m")

# This y_planned_rendezvous must be the bird's y-coordinate from t4 onward in the planned scenario.
# Let Y = y_b(t4). Since bird's planned final leg has no y-motion, Y = y_planned_rendezvous.
Y = y_planned_rendezvous
# At t5=23s, the gust hits. The bird is 2s into its planned final leg. Its y-coordinate is unchanged.
y_b_t5 = Y

# Now, we analyze the actual flight from t5 to t6.
# The bird travels from its position at t5 to the final rendezvous point (0, y_actual, 0) in (t6-t5)=17s.
# Let the bird's speed be v. The total distance squared the bird travels is D_sq = (17*v)^2.
# The distance squared can also be expressed by the change in coordinates:
# D_sq = (x_final - x_b(t5))^2 + (y_final - y_b(t5))^2 + (z_final - z_b(t5))^2
# D_sq = (0 - x_b(t5))^2 + (y_actual - y_b(t5))^2 + (0 - z_b(t5))^2
# D_sq = x_b(t5)^2 + z_b(t5)^2 + (y_actual - y_b(t5))^2

# Let's find the bird's x and z coordinates at t5.
# At t4, the bird was at (X_b(t4), Y, Z_b(t4)). To reach the y-axis in T=19s, the distance in the xz-plane is sqrt(X_b(t4)^2 + Z_b(t4)^2) = 19v.
# At t5, the bird has traveled for 2s towards the y-axis. The remaining distance in the xz-plane to the y-axis is sqrt(x_b(t5)^2 + z_b(t5)^2) = (19-2)v = 17v.
# So, x_b(t5)^2 + z_b(t5)^2 = (17v)^2.

# Substitute this into the distance squared equation:
# (17v)^2 = (17v)^2 + (y_actual - y_b(t5))^2
# This simplifies to (y_actual - y_b(t5))^2 = 0, which means y_actual = y_b(t5).
y_actual_rendezvous = y_b_t5
print(f"Actual rendezvous y-coordinate: y_actual = {y_actual_rendezvous} m")

# Step 6: Calculate the final acceleration a3
# The man's final position at t6=40s is given by:
# y_m(t6) = y_m(t5) + u_m(t5)*(t6-t5) + 0.5*a3*(t6-t5)^2
# The man must reach the actual rendezvous point.
# y_actual_rendezvous = y_m_t5 + u_m_t5 * (t6 - t5) + 0.5 * a3 * (t6 - t5)**2
# Let's solve for a3.
delta_t_final = t6 - t5
a3 = (y_actual_rendezvous - y_m_t5 - u_m_t5 * delta_t_final) / (0.5 * delta_t_final**2)

print(f"\nFinal calculation for a3:")
print(f"{y_actual_rendezvous} = {y_m_t5} + {u_m_t5} * {delta_t_final} + 0.5 * a3 * {delta_t_final}^2")
print(f"{y_actual_rendezvous} = {y_m_t5 + u_m_t5 * delta_t_final} + {0.5 * delta_t_final**2} * a3")
print(f"{y_actual_rendezvous - (y_m_t5 + u_m_t5 * delta_t_final)} = {0.5 * delta_t_final**2} * a3")
a3_numerator = y_actual_rendezvous - y_m_t5 - u_m_t5 * delta_t_final
a3_denominator = 0.5 * delta_t_final**2
print(f"a3 = {a3_numerator} / {a3_denominator}")

print(f"\nThe man's constant acceleration during this final leg of his journey is {a3} m/s^2.")
print(f"<<<{a3}>>>")
