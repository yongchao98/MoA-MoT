import math

def analyze_experiment_timescales():
    """
    Calculates and compares the key timescales in the Bi-207 decay experiment
    to determine the dominant factor for time resolution requirements.
    """
    # --- Input Parameters ---
    # The activity of the source is 1 kBq.
    activity_kBq = 1.0
    # The distance between the two detectors is 1 m, so the distance from the source to one detector is 0.5 m.
    distance_source_to_detector_m = 0.5
    # We'll use a typical energy for a conversion electron from Bi-207, which is around 1 MeV.
    electron_KE_MeV = 1.0

    # --- Constants ---
    speed_of_light_m_s = 299792458.0
    electron_rest_mass_MeV_c2 = 0.511

    # --- Step 1: Calculate the timescale set by the source activity ---
    # The activity in decays per second (Bq) is 1 kBq * 1000.
    activity_Bq = activity_kBq * 1000
    # The average time between two consecutive decays is the inverse of the activity.
    # Equation: T_decay = 1 / Activity
    avg_time_between_decays_s = 1.0 / activity_Bq

    # --- Step 2: Calculate the timescale set by the detector distance (Time of Flight) ---
    # First, find the electron's speed from its kinetic energy (KE) using the relativistic formula:
    # KE = (gamma - 1) * m*c^2
    gamma = (electron_KE_MeV / electron_rest_mass_MeV_c2) + 1
    # Then find the velocity v from gamma: gamma = 1 / sqrt(1 - v^2/c^2)
    electron_speed_m_s = speed_of_light_m_s * math.sqrt(1 - 1/gamma**2)
    # Finally, calculate the time of flight (TOF).
    # Equation: T_flight = Distance / Speed
    time_of_flight_s = distance_source_to_detector_m / electron_speed_m_s

    # --- Step 3: Print results and conclusion ---
    print("Analysis of Dominant Factor for Time Resolution")
    print("="*50)
    print("This script compares the two main timescales in the experiment.")
    print("\n1. Timescale from Source Activity:")
    print(f"   The source activity is {activity_Bq} decays/second.")
    print(f"   The average time between two separate decay events is 1 / {activity_Bq} = {avg_time_between_decays_s:.3f} seconds.")
    print(f"   This is equal to {avg_time_between_decays_s * 1e3:.1f} milliseconds.")

    print("\n2. Timescale from Detector Distance:")
    print(f"   The distance from the source to a detector is {distance_source_to_detector_m} m.")
    print(f"   The time for a {electron_KE_MeV} MeV electron to travel this distance (time of flight) is {time_of_flight_s:.2e} seconds.")
    print(f"   This is equal to {time_of_flight_s * 1e9:.2f} nanoseconds.")

    print("\n--- Conclusion ---")
    print("To individually measure electrons from separate decays, the system's time resolution must be able to distinguish events.")
    print(f"The average time between decays is ~{avg_time_between_decays_s * 1e3:.1f} ms, while the electron travel time is only ~{time_of_flight_s * 1e9:.2f} ns.")
    print("The time between decays is about a million times longer than the travel time.")
    print("The rate of decays (the activity) determines how frequently events occur and thus how likely they are to overlap in time.")
    print("A higher activity would demand a better time resolution. Therefore, the activity is the dominant factor.")
    print("="*50)

if __name__ == '__main__':
    analyze_experiment_timescales()