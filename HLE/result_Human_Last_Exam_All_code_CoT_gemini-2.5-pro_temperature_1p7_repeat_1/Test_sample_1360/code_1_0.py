# The user wants to find the dominant factor setting the minimum time resolution
# for an experiment with a 1 kBq Bi-207 source.

# Define the source activity in its base unit, Becquerels (decays per second).
# 1 kBq = 1000 Bq.
activity_Bq = 1000

# The rate of decays (activity) determines the average time interval between
# any two random, uncorrelated decay events. This is the fundamental timescale
# associated with the source activity. A shorter average time requires a better
# (smaller) time resolution to avoid random pile-up of events.
# The calculation is T = 1 / R, where R is the rate.
average_time_s = 1 / activity_Bq

# Convert seconds to milliseconds for easier interpretation.
average_time_ms = average_time_s * 1000

# Print the final equation and the result.
# This calculation demonstrates the timescale associated with the source activity (Factor D).
print(f"The source activity is {activity_Bq} decays per second (Bq).")
print("The primary time resolution requirement for a counting system is driven by the need to distinguish individual random events.")
print("We can calculate the average time between these random decays to understand the relevant timescale.")
print("\nEquation: Average Time (s) = 1 / Activity (Bq)")
print(f"Calculation: Average Time = 1 / {activity_Bq}")
print(f"Result: The average time between decays is {average_time_s} seconds, or {average_time_ms} milliseconds.")
print("\nConclusion: A higher activity would lead to a shorter average time, thus requiring a better (smaller) time resolution. Therefore, the source activity is the dominant factor that sets this requirement.")
