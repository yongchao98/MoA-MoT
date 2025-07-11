def analyze_time_resolution():
    """
    Calculates the average time between decays for a Bi 207 source
    and explains why activity is the dominant factor for time resolution.
    """
    # Activity of the source in Becquerels (Bq), where 1 Bq = 1 decay/second
    # The given activity is 1 kBq.
    activity_in_kBq = 1
    activity_in_Bq = activity_in_kBq * 1000

    # The rate of decays determines the average time between each event.
    # The minimum requirement for time resolution is to be able to distinguish
    # these consecutive, random events.
    average_time_between_decays_s = 1 / activity_in_Bq

    # Convert seconds to milliseconds for easier interpretation
    average_time_between_decays_ms = average_time_between_decays_s * 1000

    print("Analyzing the dominant factor for time resolution:")
    print("-" * 50)
    print(f"Given source activity: {activity_in_kBq} kBq")
    print(f"This corresponds to {activity_in_Bq} decays per second.")
    print("\n")
    print("The average time interval (Δt) between two consecutive decay events is the inverse of the activity rate.")
    print(f"Calculation: Δt = 1 / {activity_in_Bq} s")
    print(f"Result: Δt = {average_time_between_decays_s} s = {average_time_between_decays_ms} ms")
    print("\n")
    print("Conclusion:")
    print("To individually measure the electrons produced by these decays, the detector system's time resolution must be significantly shorter than this average interval of 1 ms to avoid 'pile-up' (mistaking two separate events for one).")
    print("While other factors like electron time-of-flight (~ns) or correlated emissions (~ns) would require even better resolution to measure, the most fundamental and dominant factor that sets the minimum requirement is the overall event rate, which is the source activity.")

analyze_time_resolution()