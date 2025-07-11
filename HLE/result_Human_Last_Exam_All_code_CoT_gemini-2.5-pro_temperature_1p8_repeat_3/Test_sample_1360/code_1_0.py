import math

# Given values
activity_kbq = 1  # in kilo-Becquerel
distance = 1.0  # in meters
speed_of_light = 3.0e8 # in m/s

# --- Step 1: Analyze the time scale set by the source activity ---
# Convert activity to Becquerel (decays per second)
activity_bq = activity_kbq * 1000

# Calculate the average time interval between two consecutive decays
# This is the reciprocal of the activity
avg_time_between_decays = 1 / activity_bq

print("Analysis of Time Scales:")
print("-----------------------")
print(f"The measured activity of the source is {activity_kbq} kBq, which is {activity_bq} decays per second.")
# The request mentions outputting each number in the final equation.
print(f"The average time (T) between decays is calculated as: T = 1 / (Activity)")
print(f"T = 1 / {activity_bq} s = {avg_time_between_decays} s")
print(f"This is equal to {avg_time_between_decays * 1000} milliseconds.")
print("\n")

# --- Step 2: Analyze the time scale set by the distance ---
# The distance to each detector from the source is half the total distance.
dist_to_detector = distance / 2
# For a rough estimate, assume electrons travel near the speed of light.
time_of_flight = dist_to_detector / speed_of_light
print(f"The time-of-flight for an electron to travel {dist_to_detector} m to a detector is approximately {time_of_flight * 1e9:.2f} nanoseconds.")
print("\n")

# --- Step 3: Conclusion ---
print("Conclusion:")
print(f"The average time between decay events is ~{avg_time_between_decays * 1000:.0f} milliseconds, while the electron time-of-flight is ~{time_of_flight * 1e9:.1f} nanoseconds.")
print(f"Since the time between decays ({avg_time_between_decays:.3f} s) is about a million times longer than the time-of-flight ({time_of_flight:.2e} s), the detector system must primarily be able to handle events separated by this average decay time.")
print("Therefore, the random rate of decays, given by the source activity, is the dominant factor setting the minimum time resolution requirement.")
