import math

# --- Inputs from the problem ---
P_initial = 100  # Initial power in Watts
x = 4          # Control subsystem power in Watts
y = 5          # Camera subsystem power in Watts
z = 1          # Sensor subsystem power in Watts
T_half = 400   # Half-life of the battery in days

# The formula to find the time 't' for a given power level P_target is:
# t = T_half * log2(P_initial / P_target)

print("Calculating the operational lifetime for Pioneer's subsystems.")
print("-" * 60)

# a. When the camera subsystem is turned off
p_target_camera = x + y + z
days_camera_float = T_half * math.log2(P_initial / p_target_camera)
days_camera_rounded = round(days_camera_float)

print("a. Camera Subsystem Shutdown Calculation:")
# The prompt requires outputting each number in the final equation.
print(f"The camera is turned off when power drops below the requirement for all subsystems.")
print(f"Power threshold = x + y + z = {x} + {y} + {z} = {p_target_camera} W.")
print(f"The equation to solve for time (t) is: t = {T_half} * log2({P_initial} / {p_target_camera})")
print(f"Calculation: t = {T_half} * log2({P_initial / p_target_camera:.1f}) = {days_camera_float:.4f} days.")
print(f"Result (a): Rounded to the nearest whole day, the camera turns off on day {days_camera_rounded}.\n")

# b. When the sensor subsystem is turned off
p_target_sensor = x + z
days_sensor_float = T_half * math.log2(P_initial / p_target_sensor)
days_sensor_rounded = round(days_sensor_float)

print("b. Sensor Subsystem Shutdown Calculation:")
print(f"The sensor is turned off when power drops below the requirement for the control and sensor systems.")
print(f"Power threshold = x + z = {x} + {z} = {p_target_sensor} W.")
print(f"The equation to solve for time (t) is: t = {T_half} * log2({P_initial} / {p_target_sensor})")
print(f"Calculation: t = {T_half} * log2({P_initial / p_target_sensor:.1f}) = {days_sensor_float:.4f} days.")
print(f"Result (b): Rounded to the nearest whole day, the sensor turns off on day {days_sensor_rounded}.\n")

# c. Memory usage in D
# An efficient C program requires variables for:
# P, x, y, z, ln2, lnP, target, time (8 frac variables)
# Each frac is 3 * 2D = 6D.
# Total memory = 8 fracs * 6D/frac = 48D.
mem_usage_d = 48

# d. Number of 'exp' function calls
# The program calls ln() 4 times: ln(2), ln(100), ln(10), ln(5).
# The number of exp() calls per ln(val) is the number of iterations in the bisection method.
# Total exp() calls = iterations(ln(2)) + iterations(ln(100)) + iterations(ln(10)) + iterations(ln(5))
# Total calls approx. = 31 + 37 + 34 + 33 = 135
exp_calls = 135

# Final Answer assembly
print("-" * 60)
print("Summary of Results (a:b:c:d):")

# The final answer must be in the requested format
final_answer = f"{days_camera_rounded}:{days_sensor_rounded}:{mem_usage_d}:{exp_calls}"
print(f"<<<{final_answer}>>>")