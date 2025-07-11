import math

# Given values
activity_in_kBq = 1  # in kilobecquerel
distance_between_detectors = 1  # in meters

# Constants
c = 299792458  # Speed of light in m/s

# 1. Calculate the decay rate (decays per second)
activity_in_Bq = activity_in_kBq * 1000

# 2. Calculate the average time between decays
# This is the most critical factor for setting the time resolution requirement.
average_time_between_decays_s = 1 / activity_in_Bq
average_time_between_decays_ms = average_time_between_decays_s * 1000
average_time_between_decays_ns = average_time_between_decays_s * 1e9

# 3. Calculate the time of flight of an electron from the source to one detector
# This is to show that it is not the dominant factor.
distance_to_detector = distance_between_detectors / 2
# Assume the electron travels near the speed of light (it's relativistic)
time_of_flight_s = distance_to_detector / c
time_of_flight_ns = time_of_flight_s * 1e9

# 4. Print the explanation
print("Analysis of Time Resolution Factors:")
print("-" * 40)
print(f"The activity of the source is {activity_in_kBq} kBq, which is equal to {activity_in_Bq} decays per second.")
print(f"This means, on average, the time between two consecutive decay events is 1 / {activity_in_Bq} = {average_time_between_decays_s:.4f} seconds, or {average_time_between_decays_ms:.1f} milliseconds.")
print("\nRadioactive decay is a random process. To distinguish one decay event from the next, the detector's time resolution must be significantly shorter than this average time.")
print(f"This average separation of {average_time_between_decays_ms:.1f} ms sets the primary scale for the required time resolution.")

print("\nLet's compare this to the time-of-flight determined by the detector distance:")
print(f"The distance from the source to a detector is {distance_to_detector} m.")
print(f"The time for a relativistic electron to travel this distance is approximately {time_of_flight_ns:.2f} nanoseconds.")

print("\nConclusion:")
print(f"The average time between decays ({average_time_between_decays_ns:,.0f} ns) is much larger than the electron's time-of-flight ({time_of_flight_ns:.2f} ns).")
print("Therefore, the dominant factor that sets the minimum time resolution requirement is the rate of decays, which is determined by the measured activity of the source.")
print("-" * 40)
print("\nThe correct option is D because the activity determines the average frequency of events that need to be individually resolved.")
