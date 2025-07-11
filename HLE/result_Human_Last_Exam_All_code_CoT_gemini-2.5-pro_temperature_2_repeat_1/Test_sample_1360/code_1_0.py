# The activity of the source is given in kilo-Becquerel (kBq).
# 1 Bq = 1 decay per second.
# 1 kBq = 1000 decays per second.

activity_in_kBq = 1
activity_in_Bq = activity_in_kBq * 1000  # Convert kBq to Bq (decays/sec)

# The time resolution requirement is determined by the average time between consecutive decay events.
# If events happen too close together, the detector cannot distinguish them.
# The average time interval is the inverse of the decay rate (activity).
average_time_interval_seconds = 1 / activity_in_Bq

# Convert the time to more readable units (milliseconds).
# 1 second = 1000 milliseconds.
average_time_interval_ms = average_time_interval_seconds * 1000

print(f"The source activity is {activity_in_Bq} decays per second.")
print(f"This means, on average, the time between two consecutive decay events is:")
print(f"1 / {activity_in_Bq} decays/second = {average_time_interval_seconds} seconds")
print(f"This is equal to {average_time_interval_ms} milliseconds.")
print("\nThis average time between decays is the dominant factor setting the time resolution requirement.")
print("The detector system must be fast enough to distinguish events that occur on this timescale.")
print("\nTherefore, the dominant factor is the measured activity of the source (1 kBq).")
