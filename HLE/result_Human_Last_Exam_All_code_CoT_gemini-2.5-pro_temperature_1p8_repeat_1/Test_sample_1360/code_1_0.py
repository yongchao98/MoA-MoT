import math

# The key to this problem is understanding the different time scales involved.
# The dominant factor is the one that sets the most critical constraint for individually measuring events.
# In a particle counting experiment, "individually measuring" primarily means avoiding event pile-up,
# which happens when two distinct events occur too close in time for the detector to resolve them.

# The rate at which events occur is given by the source activity.
activity_in_kBq = 1
activity_in_Bq = activity_in_kBq * 1000  # Convert kBq to Bq (decays per second)

# The average time between consecutive decay events is the reciprocal of the activity.
# This time scale dictates the required time resolution to keep pile-up low.
# Equation: T_avg = 1 / Activity
average_time_between_decays_s = 1 / activity_in_Bq

# For clarity, let's convert this time into milliseconds (ms) and microseconds (us).
average_time_between_decays_ms = average_time_between_decays_s * 1000
average_time_between_decays_us = average_time_between_decays_s * 1_000_000

print("Source Activity and Time Resolution Requirement Analysis:")
print("-" * 55)
print(f"The activity of the source is {activity_in_Bq} decays per second.")
print("\nThe fundamental equation for the average time (T) between decays is:")
print("T = 1 / Activity")
print(f"T = 1 / {activity_in_Bq} s")
print(f"T = {average_time_between_decays_s} s")
print("\nThis corresponds to:")
print(f"  {average_time_between_decays_ms:.1f} milliseconds (ms)")
print(f"  {average_time_between_decays_us:.1f} microseconds (µs)")
print("\nConclusion:")
print("The detector's time resolution must be significantly shorter than this average time of 1 ms")
print("to distinguish between successive decay events. Therefore, the source activity (1 kBq) is")
print("the dominant factor setting this requirement.")
