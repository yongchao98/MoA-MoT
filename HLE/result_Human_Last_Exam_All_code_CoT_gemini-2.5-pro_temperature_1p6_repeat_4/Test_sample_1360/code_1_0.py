def calculate_time_resolution_factor():
    """
    Calculates the average time between decays for a radioactive source
    and explains why it's the dominant factor for time resolution.
    """
    # Activity of the source in decays per second (1 kBq = 1000 Bq)
    activity_bq = 1000.0  # decays/second

    # The average time interval between two consecutive decays is the reciprocal of the activity.
    average_time_between_decays_s = 1.0 / activity_bq

    # Convert seconds to milliseconds for easier interpretation
    average_time_between_decays_ms = average_time_between_decays_s * 1000

    print("Source Activity: {} decays/second ({} kBq)".format(int(activity_bq), int(activity_bq/1000)))
    print("The rate of decays determines the average time between the emission of electrons.")
    print("Average time between consecutive decay events = 1 / (Activity)")
    print("Calculation: 1 / {} decays/s = {} seconds".format(int(activity_bq), average_time_between_decays_s))
    print("Average time between consecutive decay events: {} ms\n".format(average_time_between_decays_ms))

    # For comparison, let's estimate the time-of-flight
    # Distance is 0.5 m (source to one detector). Speed of electron is ~0.94c for 1 MeV electron.
    # c = 3.0e8 m/s. Let's approximate speed as c for simplicity.
    distance_m = 0.5
    speed_of_light_m_per_s = 3.0e8
    time_of_flight_s = distance_m / speed_of_light_m_per_s
    time_of_flight_ns = time_of_flight_s * 1e9 # convert to nanoseconds

    print("For comparison, the time-of-flight of an electron from the source to a detector (0.5 m away):")
    print("Time of flight is on the order of {} nanoseconds.\n".format(round(time_of_flight_ns, 2)))

    print("Conclusion:")
    print("The average time between decays (~{} ms) is much larger than the electron's time-of-flight (~{} ns).".format(
        int(average_time_between_decays_ms), round(time_of_flight_ns, 2)))
    print("To resolve individual electrons from different decays, the detector system's time resolution must be significantly better than the average time separating them.")
    print("Therefore, the dominant factor that sets the scale for the required time resolution is the rate of decays, which is the source's activity.")

calculate_time_resolution_factor()