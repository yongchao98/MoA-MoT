import math

def analyze_decay_timing():
    """
    Calculates and compares the key timescales involved in the experiment
    to determine the dominant factor for time resolution requirements.
    """
    # Given parameters
    activity_kbq = 1.0  # in kilobequerels
    distance_between_detectors_m = 1.0  # in meters
    speed_of_light_m_s = 299792458  # in meters per second

    # 1. Calculate the timescale set by the source activity
    activity_bq = activity_kbq * 1000
    average_time_between_decays_s = 1.0 / activity_bq
    average_time_between_decays_ms = average_time_between_decays_s * 1000

    # 2. Calculate the timescale set by the distance (Time-of-Flight)
    distance_to_detector_m = distance_between_detectors_m / 2.0
    # Assuming electron speed is close to the speed of light
    time_of_flight_s = distance_to_detector_m / speed_of_light_m_s
    time_of_flight_ns = time_of_flight_s * 1e9

    # 3. Print and explain the results
    print("--- Timescale Analysis ---")
    print(f"Source Activity: {int(activity_bq)} decays/second.")
    print(f"The activity determines the average time between individual decay events.")
    print(f"Average time between decays = 1 / {int(activity_bq)} Hz = {average_time_between_decays_s:.3f} seconds = {average_time_between_decays_ms:.1f} ms.")
    print("\n")
    print(f"The distance to each detector is {distance_to_detector_m} m.")
    print("This determines the electron's time-of-flight (TOF).")
    print(f"Time-of-Flight ≈ {distance_to_detector_m} m / c ≈ {time_of_flight_ns:.2f} nanoseconds.")
    print("\n--- Conclusion ---")
    print(f"The time between decay events (~{average_time_between_decays_ms:.1f} ms) is the dominant timescale.")
    print("This is the timescale the system must resolve to distinguish one independent event from the next.")
    print("It is approximately a million times longer than the time-of-flight.")
    print("Therefore, the dominant factor setting the minimum requirement to measure events individually is the source activity.")

analyze_decay_timing()
