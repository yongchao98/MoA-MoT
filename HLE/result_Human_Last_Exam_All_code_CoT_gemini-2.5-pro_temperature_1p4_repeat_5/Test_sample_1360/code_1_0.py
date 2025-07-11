import math

# --- Given information ---
# Choice D: Activity is given in kilo-Becquerel (kBq)
activity_in_kBq = 1
# Choice A: Distance between detectors in meters
distance_detectors = 1.0

# --- Analysis ---

# 1. Analyze the effect of Activity (Choice D)
# Activity is the number of decays per second.
# 1 kBq = 1000 Bq (decays per second)
activity_in_Bq = activity_in_kBq * 1000

# The average time interval between two random, independent decay events is the
# reciprocal of the activity rate. To resolve two separate decay events,
# the detector's time resolution must be smaller than this average interval.
# The calculation is: Average Time = 1 / Activity
numerator_for_time = 1
denominator_for_time = activity_in_Bq
average_time_between_decays_s = numerator_for_time / denominator_for_time

print("--- Analysis of Factor: Source Activity (Choice D) ---")
print(f"The source activity is {activity_in_Bq} decays per second.")
print(f"The average time between decays is calculated by the equation: {numerator_for_time} / {denominator_for_time}")
print(f"Result: The average time between decay events is {average_time_between_decays_s:.3f} seconds, or {average_time_between_decays_s * 1000:.1f} milliseconds.\n")


# 2. Analyze the effect of Distance (Choice A)
# This determines the time-of-flight (TOF) of an electron from source to detector.
# The source is in the middle, so the distance is half the total.
distance_to_detector_m = distance_detectors / 2

# We need the electron's speed. Bi-207 emits electrons up to ~1 MeV.
# We can calculate the speed of a 1 MeV electron, which is close to the speed of light.
# Physical constants
c = 299792458  # Speed of light in m/s
m0c2_electron_MeV = 0.511  # Electron rest mass energy in MeV

# Electron kinetic energy in MeV
K_electron_MeV = 1.0
# Total energy E = gamma * m0c2
E_total_MeV = K_electron_MeV + m0c2_electron_MeV
# Lorentz factor gamma
gamma = E_total_MeV / m0c2_electron_MeV
# Electron speed v = c * sqrt(1 - 1/gamma^2)
v_electron_ms = c * math.sqrt(1 - 1 / (gamma**2))
# Time of flight = distance / speed
time_of_flight_s = distance_to_detector_m / v_electron_ms

print("--- Analysis of Factor: Distance (Choice A) ---")
print(f"The distance from source to one detector is {distance_to_detector_m} meters.")
print(f"The calculated time-of-flight for a high-energy electron is {time_of_flight_s * 1e9:.2f} nanoseconds.\n")

# 3. Conclusion
print("--- Conclusion ---")
print("To individually measure electrons from separate decays, we must distinguish one decay event from the next.")
print(f"The timescale determined by activity is ~{average_time_between_decays_s * 1000:.1f} milliseconds.")
print(f"The timescale determined by distance is ~{time_of_flight_s * 1e9:.2f} nanoseconds.")
print(f"The average time between decays ({average_time_between_decays_s * 1000:.1f} ms) is hundreds of thousands of times longer than the time of flight ({time_of_flight_s * 1e9:.2f} ns).")
print("The rate of decays (the activity) dictates how often events occur. This sets the primary requirement for time resolution to prevent signals from different decays from overlapping (pile-up). Therefore, it is the dominant factor.")
