import math

# Plan:
# 1. Define physical constants from the problem.
# 2. Calculate the moth's speed relative to the ground.
# 3. Calculate the time until the LED blink sequence is triggered.
# 4. Calculate the time it takes for the blink sequence to reach the last LED.
# 5. Sum these times to get the total time elapsed from the moth's start.
# 6. Calculate the moth's displacement over this total time.
# 7. Print the calculation steps and the final answer, showing the numbers in the equation.

# 1. Define constants
air_speed_m_per_min = 5.0
moth_airspeed_m_per_min = 5.675
num_leds = 80
led_blink_delay_s = 0.3
moth_start_pos_m = 2.0
trigger_pos_m = 1.0

print("Step 1: Calculate the moth's velocity relative to the ground.")
# The moth flies west (negative direction) and the wind blows east (positive direction).
v_moth_ground_m_per_min = -moth_airspeed_m_per_min + air_speed_m_per_min
print(f"Moth's ground velocity = (Air Speed) - (Moth's Airspeed) = {air_speed_m_per_min} m/min - {moth_airspeed_m_per_min} m/min = {v_moth_ground_m_per_min:.3f} m/min.")

print("\nStep 2: Calculate the time until the LED sequence starts.")
distance_to_trigger_m = moth_start_pos_m - trigger_pos_m
speed_moth_ground_m_per_min = abs(v_moth_ground_m_per_min)
time_to_trigger_min = distance_to_trigger_m / speed_moth_ground_m_per_min
print(f"The moth must travel {distance_to_trigger_m} m at a speed of {speed_moth_ground_m_per_min:.3f} m/min.")
print(f"Time to trigger = Distance / Speed = {distance_to_trigger_m} m / {speed_moth_ground_m_per_min:.3f} m/min = {time_to_trigger_min:.4f} min.")

print("\nStep 3: Calculate the time for the light signal to travel from the first to the last LED.")
intervals = num_leds - 1
cascade_time_s = intervals * led_blink_delay_s
cascade_time_min = cascade_time_s / 60.0
print(f"The signal travels across {intervals} intervals with a {led_blink_delay_s} s delay each.")
print(f"Cascade time = {intervals} * {led_blink_delay_s} s = {cascade_time_s:.1f} s, which is {cascade_time_min:.4f} min.")

print("\nStep 4: Calculate the final displacement.")
total_time_min = time_to_trigger_min + cascade_time_min
total_displacement_m = v_moth_ground_m_per_min * total_time_min
print(f"The total time elapsed from the start is the sum of the time to trigger and the cascade time.")
print(f"Total time = {time_to_trigger_min:.4f} min + {cascade_time_min:.4f} min = {total_time_min:.4f} min.")
print("\nThe final displacement is the moth's ground velocity multiplied by the total time.")
print(f"Displacement = (Velocity) * (Total Time)")
print(f"Displacement = ({v_moth_ground_m_per_min:.3f} m/min) * ({total_time_min:.4f} min)")
print(f"Displacement = {total_displacement_m:.4f} m")

print(f"\nThe moth's displacement is a vector of {total_displacement_m:.4f} m. The magnitude (distance) of this displacement is {abs(total_displacement_m):.2f} m.")
