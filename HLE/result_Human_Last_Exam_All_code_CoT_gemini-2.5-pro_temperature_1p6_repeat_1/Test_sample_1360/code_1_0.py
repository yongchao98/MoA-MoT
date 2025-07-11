import math

def analyze_time_resolution():
    """
    Calculates and compares timescales in the Bi-207 decay experiment
    to determine the dominant factor for time resolution.
    """
    # Given values
    activity_bq = 1000  # 1 kBq = 1000 decays/second
    distance_to_detector_m = 0.5  # Source is in the middle, 1 m apart
    speed_of_light_mps = 3.0e8

    # --- Step 1: Calculate the timescale set by the source activity ---
    # The activity determines the average rate of events. The average time
    # between these events is the reciprocal of the activity.
    avg_time_between_decays_s = 1.0 / activity_bq

    print("--- Analysis of Timescales ---")
    print(f"The source activity is {activity_bq} decays per second.")
    print("This means that, on average, a decay event occurs every:")
    print(f"  T_avg = 1 / {activity_bq} Hz = {avg_time_between_decays_s} seconds")
    print(f"  T_avg = {avg_time_between_decays_s * 1000:.1f} milliseconds (ms)")
    print(f"  T_avg = {avg_time_between_decays_s * 1e6:.1f} microseconds (µs)")
    print("\nThis average time between decays is the characteristic time of the event stream. To measure events individually and avoid pile-up, the detector system's time resolution must be significantly shorter than this duration.")
    print("-" * 35)

    # --- Step 2: Calculate another relevant timescale (electron time-of-flight) ---
    # The energy of conversion electrons from Bi-207 is up to ~1 MeV.
    # An electron with 1 MeV kinetic energy travels at about 94% of the speed of light.
    electron_velocity_mps = 0.94 * speed_of_light_mps
    time_of_flight_s = distance_to_detector_m / electron_velocity_mps

    print("Now, let's consider the timescale set by the distance between detectors.")
    print("The time for an electron to travel 0.5 m to a detector (Time-of-Flight) is approximately:")
    print(f"  t_flight = 0.5 m / (0.94 * c) = {time_of_flight_s:.3e} seconds")
    print(f"  t_flight = {time_of_flight_s * 1e9:.2f} nanoseconds (ns)")
    print("\nOther phenomena, like correlated emissions from a single decay, also occur on nanosecond or shorter timescales.")
    print("-" * 35)

    # --- Step 3: Compare the timescales and conclude ---
    print("--- Conclusion ---")
    print(f"Comparing the timescales:")
    print(f"  - Time between decays (from Activity): {avg_time_between_decays_s * 1e6:,.0f} µs")
    print(f"  - Electron time-of-flight (from Distance): {time_of_flight_s * 1e6:.4f} µs")
    print("\nThe time between separate decay events is orders of magnitude longer than the electron's travel time.")
    print("Therefore, the dominant factor that sets the minimum requirement for telling one event from the next is the rate at which they occur, which is the source's activity.")
    print("A higher activity would lead to a shorter average time between events, thus requiring a better (faster) time resolution.")

if __name__ == '__main__':
    analyze_time_resolution()
<<<D>>>