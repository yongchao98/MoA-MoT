import math

# Man's initial parameters
u0 = 1.5  # m/s
t0 = 0
t1 = 4    # s

# Segment 1: t0 to t1 (Constant speed)
y_t1 = u0 * t1
u_t1 = u0

# Segment 2: t1 to t2 (Deceleration)
a1 = -0.15 # m/s^2
# Time to stop (u=0): u_t1 + a1*(t2-t1) = 0 => 1.5 - 0.15*(t2-4) = 0
# t2-4 = 1.5/0.15 = 10 => t2 = 14s
t2 = 14
dt_12 = t2 - t1
y_t2 = y_t1 + u_t1 * dt_12 + 0.5 * a1 * dt_12**2
u_t2 = 0

# Segment 3: t2 to t3 (Still)
t3 = 15   # s
y_t3 = y_t2
u_t3 = 0

# Segment 4: t3 to t4 (Acceleration)
a2 = 0.25 # m/s^2
# Speed at t4 is u0=1.5 m/s
# u_t4 = u_t3 + a2*(t4-t3) => 1.5 = 0 + 0.25*(t4-15)
# t4-15 = 1.5/0.25 = 6 => t4 = 21s
t4 = 21
dt_34 = t4 - t3
y_t4 = y_t3 + u_t3 * dt_34 + 0.5 * a2 * dt_34**2
u_t4 = u_t3 + a2 * dt_34

# Segment 5 (original plan): t4 onwards (Constant speed u_t4)
# In the original plan, the man continues at this constant speed.
# Let's find the y position at t_orig = 40s in this original plan.
t_orig = 40
dt_4_orig = t_orig - t4
y_orig_40 = y_t4 + u_t4 * dt_4_orig

# This y_orig_40 is the y-coordinate of the new rendezvous point.
y_new = y_orig_40

# Now, let's analyze the man's actual motion from t5 onwards to find a3.
# First, find man's state at t5.
t5 = 23   # s
dt_45 = t5 - t4
y_t5 = y_t4 + u_t4 * dt_45
u_t5 = u_t4

# Segment 6 (actual plan): t5 to t6 (Acceleration a3)
t6 = 40   # s
dt_56 = t6 - t5
# Final position y_new is given by: y_new = y_t5 + u_t5*dt_56 + 0.5*a3*dt_56**2
# We can solve for a3:
# a3 = (y_new - y_t5 - u_t5*dt_56) * 2 / dt_56**2
a3 = (y_new - y_t5 - u_t5 * dt_56) * 2 / (dt_56**2)

# Print the calculation steps
print(f"Man's position at t1 = {t1} s is y_t1 = {y_t1:.2f} m.")
print(f"Man stops at t2 = {t2} s. Position is y_t2 = {y_t2:.2f} m.")
print(f"Man is at rest until t3 = {t3} s. Position is y_t3 = {y_t3:.2f} m.")
print(f"Man accelerates until t4 = {t4} s. Position is y_t4 = {y_t4:.2f} m, speed is u_t4 = {u_t4:.2f} m/s.")
print(f"Original rendezvous at y_orig = {y_orig_40:.2f} m.")
print(f"New rendezvous is at y_new = {y_new:.2f} m.")
print(f"Man's state at t5 = {t5} s: position y_t5 = {y_t5:.2f} m, speed u_t5 = {u_t5:.2f} m/s.")
print(f"The final leg is from t5 to t6 = {t6} s, a duration of {dt_56} s.")

print(f"Using the equation for motion: y_new = y_t5 + u_t5*dt_56 + 0.5*a3*dt_56^2")
print(f"{y_new:.2f} = {y_t5:.2f} + {u_t5:.2f} * {dt_56} + 0.5 * a3 * {dt_56**2}")
print(f"{y_new:.2f} = {y_t5 + u_t5 * dt_56:.2f} + {0.5 * dt_56**2} * a3")
print(f"{y_new - (y_t5 + u_t5 * dt_56):.2f} = {0.5 * dt_56**2} * a3")
print(f"a3 = {(y_new - (y_t5 + u_t5 * dt_56)) * 2} / {dt_56**2}")
print(f"a3 = {a3:.4f} m/s^2")

# Final rounded answer
final_a3 = round(a3, 2)
print(f"\nThe man's constant acceleration during this final leg is {final_a3} m/s^2.")
print(f"<<<{final_a3}>>>")