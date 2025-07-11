def calculate_time_resolution_factor():
    """
    Calculates the average time between decays for a radioactive source
    and explains why this is the dominant factor for time resolution.
    """
    # Activity in decays per second (1 kBq = 1000 Bq)
    activity_bq = 1000

    # Calculate the average time interval between decays
    average_time_interval_s = 1 / activity_bq
    average_time_interval_ms = average_time_interval_s * 1000

    print("The activity of the source is 1 kBq, which is 1000 decays per second.")
    print("To measure electrons from individual decay events, we must first determine the average time between these events.")
    print("The calculation for the average time interval (Δt) is: Δt = 1 / Activity")
    print("\n--- Calculation ---")
    print(f"Δt = 1 / {activity_bq} decays/s = {average_time_interval_s} s")
    print(f"This is equal to {average_time_interval_ms} milliseconds.")
    print("\n--- Conclusion ---")
    print("The detector system's time resolution must be significantly better than this 1 ms interval to distinguish one decay event from the next.")
    print("Other factors, like the electron's time-of-flight (~1.7 nanoseconds), are millions of times shorter and thus not the dominant constraint.")
    print("Therefore, the dominant factor that sets the minimum time resolution requirement is the measured activity of the source.")

calculate_time_resolution_factor()