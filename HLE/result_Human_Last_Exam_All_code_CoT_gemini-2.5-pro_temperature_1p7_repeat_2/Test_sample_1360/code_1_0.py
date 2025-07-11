# The activity of the source in Becquerel (Bq), which is decays per second.
# 1 kBq = 1000 Bq.
activity_Bq = 1000

# The average time between two consecutive random decays is the inverse of the activity.
# This determines the timescale requirement for the detector system.
average_time_s = 1 / activity_Bq

# Convert the time to milliseconds for easier interpretation.
average_time_ms = average_time_s * 1000

print(f"The activity of the source is {activity_Bq} decays/second.")
print(f"The average time between decays is calculated as: 1 / {activity_Bq} decays/second = {average_time_s} seconds.")
print(f"This is equal to {average_time_ms} milliseconds.")
print("\nTo measure individual decay events, the detector's time resolution must be much shorter than this average time separation.")
print("Therefore, the measured activity of the source is the dominant factor setting this requirement.")
