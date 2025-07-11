# The activity of the source in kiloBecquerels (kBq)
activity_kBq = 1

# Convert activity to Becquerels (decays per second)
# 1 kBq = 1000 Bq
activity_Bq = activity_kBq * 1000

# The average time between decays is the reciprocal of the decay rate (activity)
average_time_seconds = 1 / activity_Bq

# Convert the time to milliseconds for easier interpretation
average_time_ms = average_time_seconds * 1000

print(f"Source activity: {activity_kBq} kBq")
print(f"This corresponds to {activity_Bq} decays per second.")
print(f"The average time between two consecutive decay events is {average_time_seconds} seconds, or {average_time_ms} milliseconds.")
print("\nConclusion:")
print("To measure electrons from individual decays, the detector system's time resolution must be able to distinguish events that are, on average, 1 ms apart.")
print("This time scale (ms) is determined by the source activity and is much longer than other time scales in the experiment (like the electron's nanosecond time-of-flight).")
print("Therefore, the dominant factor that sets the time resolution requirement is the source activity.")
