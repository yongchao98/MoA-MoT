# The activity of the source is given in kilo-Becquerel (kBq).
# 1 Bq = 1 decay per second.
# 1 kBq = 1000 Bq.
activity_in_kBq = 1
activity_in_Bq = activity_in_kBq * 1000  # Activity in decays per second

# The average time interval between random, independent events (like radioactive decays)
# is the reciprocal of the event rate (the activity).
# Formula: Average Time = 1 / Activity
average_time_seconds = 1 / activity_in_Bq

# Convert seconds to milliseconds for easier interpretation.
# 1 second = 1000 milliseconds.
average_time_ms = average_time_seconds * 1000

print("Given the source activity, we can calculate the average time between consecutive decay events.")
print(f"The equation is: Average Time = 1 / Activity")
print(f"Using the given values: Average Time = 1 / {activity_in_Bq} decays/sec")
print(f"The average time between decays is {average_time_seconds} seconds, or {average_time_ms} milliseconds.")
print("\nThis means that to distinguish between the random decay events, the detector system's time resolution must be significantly shorter than this average time. Therefore, the source activity is the dominant factor setting this requirement.")
