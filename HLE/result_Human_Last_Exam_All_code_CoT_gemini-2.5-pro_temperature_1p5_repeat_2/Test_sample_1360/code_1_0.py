import math

def calculate_time_resolution_factors():
    """
    Calculates and compares the time scales relevant to the experiment
    to determine the dominant factor for time resolution.
    """

    # --- Factor 1: Measured Activity of the Source ---
    # The source has a measured activity of 1 kBq.
    # 1 kBq = 1000 Bq = 1000 decays per second.
    activity_in_Bq = 1000

    # The average time interval between independent decay events is the reciprocal of the activity.
    # Equation: T_avg = 1 / Activity
    avg_time_between_decays_sec = 1 / activity_in_Bq

    # Convert to milliseconds for easier interpretation.
    avg_time_between_decays_ms = avg_time_between_decays_sec * 1000

    print("--- Analysis based on Source Activity ---")
    print(f"The source activity is {activity_in_Bq} decays/second.")
    print("The average time between individual decay events is calculated as: 1 / {activity_in_Bq} s".format(activity_in_Bq=activity_in_Bq))
    print(f"Result: {avg_time_between_decays_sec:.3f} seconds, or {avg_time_between_decays_ms:.1f} milliseconds.\n")

    # --- Factor 2: Distance between detectors ---
    # The distance to one detector is half the total distance.
    distance_to_detector_m = 1.0 / 2.0

    # Speed of light in m/s
    c = 299792458

    # The emitted electrons have kinetic energies up to ~1 MeV.
    # For a ~1 MeV electron, the velocity is ~94% of the speed of light.
    # v ≈ 0.94 * c
    electron_velocity = 0.94 * c

    # Calculate the Time of Flight (ToF).
    # Equation: ToF = Distance / Velocity
    time_of_flight_sec = distance_to_detector_m / electron_velocity

    # Convert to nanoseconds.
    time_of_flight_ns = time_of_flight_sec * 1e9

    print("--- Analysis based on Detector Distance ---")
    print(f"The distance from the source to a detector is {distance_to_detector_m} meters.")
    print("The approximate electron Time of Flight (ToF) is calculated as: {dist:.1f} m / ({vel_frac:.2f} * c) m/s".format(dist=distance_to_detector_m, vel_frac=0.94))
    print(f"Result: Approximately {time_of_flight_ns:.2f} nanoseconds.\n")

    # --- Conclusion ---
    print("--- Conclusion ---")
    print(f"Average time between decays: {avg_time_between_decays_ms:.1f} ms")
    print(f"Electron time of flight:    {time_of_flight_ns:.2f} ns")
    print("\nTo individually measure electrons from separate decays, the system's time resolution must be significantly shorter than the average time between them (1 ms).")
    print("The time between decays (1,000,000 ns) is much longer than the time of flight (~1.8 ns).")
    print("Therefore, the dominant factor setting this requirement is the rate at which decays occur, which is defined by the source's activity.")

if __name__ == "__main__":
    calculate_time_resolution_factors()
<<<D>>>