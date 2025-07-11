import math

# --- Step 1: Analyze the time scale from the source activity ---

# The activity is given as 1 kBq (kilo-Becquerel).
activity_in_kBq = 1

# Convert activity to Bq (Becquerels, or decays per second).
# 1 kBq = 1000 Bq.
activity_in_Bq = activity_in_kBq * 1000

# Radioactive decay is a random process. The activity represents the average rate of decay.
# The average time interval (Δt) between two consecutive, random events is the reciprocal of the event rate.
average_time_between_decays_s = 1 / activity_in_Bq

# Convert the result to a more intuitive unit, like milliseconds (ms).
# 1 s = 1000 ms.
average_time_between_decays_ms = average_time_between_decays_s * 1000

print("Step 1: Calculate the average time between decay events.")
print(f"The source activity is {activity_in_kBq} kBq, which is equal to {activity_in_Bq} decays per second.")
print("The average time interval (Δt) between decays is calculated as: 1 / Activity")
print(f"Δt = 1 / {activity_in_Bq} s = {average_time_between_decays_s} s")
print(f"This is equal to {average_time_between_decays_ms} milliseconds.\n")


# --- Step 2: Analyze the time scale from the distance ---

# The detectors are 1 m apart, with the source in the middle.
# So, the distance from the source to each detector is 0.5 m.
distance_m = 0.5

# Electrons emitted from Bi-207 decay are energetic. Their speed will be close to the speed of light.
# Speed of light (c) is approximately 3.0 x 10^8 m/s.
speed_of_light_m_per_s = 3.0e8

# Calculate the time-of-flight (ToF) for an electron.
tof_s = distance_m / speed_of_light_m_per_s

# Convert ToF to nanoseconds (ns) for easier comparison.
# 1 s = 1e9 ns.
tof_ns = tof_s * 1e9

print("Step 2: Calculate the electron time-of-flight (for comparison).")
print(f"The distance from the source to a detector is {distance_m} m.")
print(f"The time-of-flight (ToF) for a relativistic electron is: Distance / Speed of Light")
print(f"ToF = {distance_m} m / {speed_of_light_m_per_s:.1e} m/s = {tof_s:.2e} s")
print(f"This is equal to {tof_ns:.2f} nanoseconds.\n")


# --- Step 3: Compare and conclude ---

print("Conclusion:")
print(f"The average time between decays is {average_time_between_decays_ms} ms, while the electron travel time is only about {tof_ns:.2f} ns.")
print(f"Since 1 millisecond is 1,000,000 nanoseconds, the time interval between decays is much larger than the travel time.")
print("To measure electrons individually, the system must be able to resolve events that are, on average, 1 ms apart.")
print("Therefore, the rate of decays (the activity) is the dominant factor setting the minimum time resolution requirement.")