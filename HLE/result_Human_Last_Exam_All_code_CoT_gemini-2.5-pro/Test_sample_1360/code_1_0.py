import math

def solve_time_resolution():
    """
    Calculates the dominant factor for time resolution in a nuclear decay experiment.
    """
    # Given parameters
    activity_kBq = 1.0  # in kilo-Becquerel
    distance_m = 1.0  # distance between detectors

    # Convert activity to decays per second (Becquerel, Bq)
    activity_Bq = activity_kBq * 1000

    # --- Step 1: Calculate the average time between decays ---
    # The activity (A) is the number of decays per unit time (N/t).
    # The average time interval (Δt) between individual decays is the inverse of the activity.
    # This value dictates the time resolution needed to distinguish one decay from the next.
    avg_time_between_decays_s = 1.0 / activity_Bq

    print("The dominant factor is the source activity because it determines the rate of events.")
    print("Let's calculate the average time between decay events based on the activity.")
    print("-" * 60)
    print(f"The source activity is {activity_kBq} kBq, which is {activity_Bq:.0f} decays per second.")
    print("\nThe average time interval (Δt) between decays is calculated as:")
    print(f"Δt = 1 / Activity")
    # Final equation with numbers as requested
    print(f"Δt = 1 / {activity_Bq:.0f} decays/second = {avg_time_between_decays_s} seconds")
    print(f"This is equal to {avg_time_between_decays_s * 1000:.1f} milliseconds (ms).\n")
    print("A detector system must have a time resolution significantly better than this value")
    print("to be able to resolve most decay events as individual signals.")
    print("-" * 60)

    # --- Step 2 (for comparison): Calculate the electron's time of flight ---
    # This shows why distance is not the dominant factor for resolving *separate* decays.
    # The kinetic energy of Bi-207 conversion electrons is ~1 MeV.
    electron_KE_MeV = 1.0
    electron_rest_mass_MeV = 0.511
    speed_of_light_m_s = 299792458

    # Relativistic calculation for velocity
    gamma = (electron_KE_MeV / electron_rest_mass_MeV) + 1
    beta = math.sqrt(1 - 1 / (gamma**2))
    velocity_m_s = beta * speed_of_light_m_s
    distance_to_detector_m = distance_m / 2.0
    time_of_flight_s = distance_to_detector_m / velocity_m_s

    print("For comparison, the time for an electron to travel to a detector is:")
    print(f"Time of Flight = Distance / Velocity = {distance_to_detector_m:.1f} m / {velocity_m_s:.0f} m/s = {time_of_flight_s:.3e} seconds")
    print(f"This is equal to {time_of_flight_s * 1e9:.2f} nanoseconds (ns).\n")
    print(f"Conclusion: The time between decays ({avg_time_between_decays_s * 1000:.1f} ms) is about a million times larger")
    print(f"than the electron travel time ({time_of_flight_s * 1e9:.2f} ns). Therefore, the activity is the dominant factor.")

solve_time_resolution()
<<<D>>>