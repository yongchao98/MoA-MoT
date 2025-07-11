def calculate_time_resolution_factor():
    """
    Calculates the average time between radioactive decays for a given activity
    and explains why activity is the dominant factor for time resolution.
    """
    # Activity of the Bi 207 source in Becquerels (Bq)
    # 1 kBq = 1000 Bq (decays per second)
    activity_bq = 1000

    # The average time interval between decays is the inverse of the activity.
    # This sets the timescale for resolving individual decay events.
    average_time_between_decays_s = 1 / activity_bq

    # Convert the time to milliseconds for better readability
    average_time_between_decays_ms = average_time_between_decays_s * 1000

    print("Given the source activity of 1 kBq (1000 decays/second):")
    print(f"The average time between individual decay events is 1 / {activity_bq} s = {average_time_between_decays_s} seconds.")
    print(f"This is equal to {average_time_between_decays_ms} milliseconds.")
    print("\nTo measure electrons from each decay individually, the detector system's time resolution")
    print("must be significantly shorter than this average time interval.")
    print("Other factors, like the electron's time-of-flight (~1.7 nanoseconds), are much smaller")
    print("and therefore not the dominant constraint.")
    print("\nConclusion: The dominant factor that sets the minimum time resolution requirement is the")
    print("rate of decays, which is the measured activity of the source.")

calculate_time_resolution_factor()