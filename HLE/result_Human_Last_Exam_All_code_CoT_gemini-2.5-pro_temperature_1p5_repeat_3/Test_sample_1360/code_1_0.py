import math

def solve_time_resolution():
    """
    This function calculates and compares the key timescales in the experiment
    to determine the dominant factor for time resolution.
    """

    # 1. Analyze the timescale dictated by the source activity.
    activity_in_Bq = 1000  # 1 kBq = 1000 decays/second

    # The average time interval between decays is the inverse of the activity.
    avg_time_between_decays_s = 1 / activity_in_Bq
    avg_time_between_decays_ms = avg_time_between_decays_s * 1000

    print("Step 1: Calculate the average time between consecutive decay events.")
    print("The source activity determines the rate at which decay events occur.")
    print(f"Source Activity = {activity_in_Bq} decays/second")
    print("Average time between decays = 1 / (Activity)")
    print(f"                                = 1 / {activity_in_Bq} s")
    print(f"                                = {avg_time_between_decays_s} seconds, or {avg_time_between_decays_ms:.1f} milliseconds (ms).\n")

    # 2. Analyze the timescale dictated by the distance between detectors.
    distance_to_detector_m = 0.5  # Source is in the middle of 1m span
    speed_of_light_mps = 299792458  # m/s

    # The time of flight is distance / speed.
    # Electrons from Bi-207 decay are relativistic, with speeds close to c.
    time_of_flight_s = distance_to_detector_m / speed_of_light_mps
    time_of_flight_ns = time_of_flight_s * 1e9

    print("Step 2: Calculate the electron's time of flight from source to detector.")
    print("This is determined by the distance and the electron's speed.")
    print(f"Distance to detector = {distance_to_detector_m} m")
    print(f"Electron speed (approx.) = Speed of light ≈ {speed_of_light_mps:.3e} m/s")
    print("Time of flight = Distance / Speed")
    print(f"               = {distance_to_detector_m} m / {speed_of_light_mps:.3e} m/s")
    print(f"               = {time_of_flight_s:.3e} seconds, or {time_of_flight_ns:.2f} nanoseconds (ns).\n")
    
    # 3. Compare the timescales and conclude.
    print("Step 3: Compare the timescales.")
    print(f"Average time between decay events: {avg_time_between_decays_ms:.1f} ms")
    print(f"Electron time of flight:             {time_of_flight_ns:.2f} ns")
    
    ratio = avg_time_between_decays_s / time_of_flight_s
    print(f"\nThe time between decays is ~{ratio:,.0f} times longer than the time of flight.")
    print("\nTo measure electrons individually (i.e., resolve one decay from the next), the detector's resolving time must be significantly shorter than the average time between the events themselves.")
    print(f"Therefore, the minimum time resolution requirement is set by the ~{avg_time_between_decays_ms:.1f} ms average separation between decays, which is a direct consequence of the source's activity.")

solve_time_resolution()
<<<D>>>