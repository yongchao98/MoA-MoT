# The activity of the source determines the average rate of decays.
# This rate is the dominant factor in setting the time resolution required
# to distinguish one decay event from the next.

# Activity in kiloBecquerel (kBq)
activity_kBq = 1

# Convert activity to Becquerel (decays per second)
# 1 kBq = 1000 Bq
activity_Bq = activity_kBq * 1000

# The average time interval (T) between consecutive decays is the reciprocal of the activity (A).
# Formula: T = 1 / A
average_time_s = 1 / activity_Bq

# Convert seconds to milliseconds for easier interpretation
average_time_ms = average_time_s * 1000

print("The activity of the source is given as 1 kBq, which is 1000 decays per second.")
print("The average time interval 'T' between two consecutive decay events is calculated as the reciprocal of the activity 'A'.")
print("\nFinal Equation: T = 1 / A")
print(f"Plugging in the numbers: T = 1 / {activity_Bq} decays/second")
print(f"\nThe calculated average time between decays is {average_time_s} seconds, or {average_time_ms} milliseconds.")
print("\nTo resolve individual decay events, the detector system's time resolution must be significantly shorter than this 1 ms timescale.")
print("This is the dominant limiting factor compared to other effects like electron time-of-flight (~1.7 ns).")
print("Therefore, the measured activity of the source is the dominant factor.")
