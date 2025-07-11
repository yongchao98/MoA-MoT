# The activity of the source is given in kilo-Becquerels (kBq).
activity_kbq = 1

# Convert activity from kBq to Bq (decays per second).
# 1 kBq = 1000 Bq.
activity_bq = activity_kbq * 1000

# The activity represents the average number of decays per second.
# The average time interval between two consecutive decays is the reciprocal of the activity.
time_interval_s = 1 / activity_bq

# Convert the time interval to more convenient units (milliseconds and microseconds).
time_interval_ms = time_interval_s * 1000
time_interval_us = time_interval_s * 1_000_000

print(f"Source Activity: {activity_kbq} kBq")
print(f"This is equal to {activity_bq} decays per second.")
print("\nThe minimum time resolution requirement is set by the need to distinguish one random decay event from the next.")
print("We can calculate the average time between these random events.")
print("\nFinal Equation:")
# The problem requests to output each number in the final equation.
print(f"Average time interval = 1 / {activity_bq} decays/s = {time_interval_s} seconds")

print(f"\nThis average time is {time_interval_ms:.1f} milliseconds or {time_interval_us:.0f} microseconds.")
print("This time scale is much larger than any other time scale in the problem (like the ~1.7 nanosecond electron time-of-flight).")
print("Therefore, the activity, which sets this average time between events, is the dominant factor.")
