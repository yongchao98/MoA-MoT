import math

def calculate_time_resolution_basis():
    """
    Calculates the average time between decays for the given source activity
    and explains why activity is the dominant factor for time resolution.
    """
    # Given activity in kilobecquerels (kBq)
    activity_kbq = 1

    # Convert activity to Becquerels (decays per second)
    activity_bq = activity_kbq * 1000

    # The rate of decay sets the average time between individual decay events.
    # To measure events individually, the time resolution of the detector system
    # must be significantly shorter than this average time separation to avoid pile-up.
    average_time_separation_s = 1 / activity_bq

    # Convert the time to milliseconds for easier interpretation
    average_time_separation_ms = average_time_separation_s * 1000

    print("The primary goal is to measure individual electrons from radioactive decays.")
    print("Radioactive decay is a random process. The key factor limiting our ability to distinguish one decay event from the next is the average time between them.")
    print("\nThis average time is determined by the source's activity.")
    print(f"The given activity is {activity_kbq} kBq, which is {activity_bq} decays per second.")
    print("\nThe calculation for the average time separation is:")
    print("Average Time = 1 / Activity")
    print(f"Average Time = 1 / {activity_bq} decays/second")
    # Using math.trunc to show integer value for the final ms result
    print(f"Average Time = {average_time_separation_s:.3f} seconds = {math.trunc(average_time_separation_ms)} milliseconds")
    print("\nTherefore, the time resolution of the measurement system must be much less than 1 ms. This requirement is set by the activity of the source (1 kBq).")

calculate_time_resolution_basis()