import math

def analyze_time_resolution():
    """
    This script analyzes the factors affecting time resolution for measuring
    electron emissions from a Bi-207 source.
    """

    # The activity of the source is given in kiloBecquerel (kBq).
    activity_kbq = 1.0

    # Convert activity to Becquerels (Bq), which are decays per second.
    activity_bq = activity_kbq * 1000

    # The rate of decays determines how often, on average, an electron is emitted.
    # The average time interval between two successive decays is the reciprocal of the activity.
    avg_time_between_decays_s = 1.0 / activity_bq

    # --- Main Analysis ---
    print("Analyzing the dominant factor for time resolution requirements.")
    print("-" * 60)
    print("The goal is to 'individually measure' electrons, which means we need to distinguish electrons from different decay events.")

    # 1. Analyze the effect of Source Activity
    print("\nStep 1: Calculate the average time between decay events from the source activity.")
    print(f"The source activity is {activity_kbq} kBq, which is {activity_bq} decays per second.")
    print("The average time between decays is calculated as: 1 / (Activity in Bq)")
    # The final equation with numbers
    print(f"Average Time = 1 / {int(activity_bq)} = {avg_time_between_decays_s} seconds")
    print(f"This is equivalent to {avg_time_between_decays_s * 1000} milliseconds (ms).")
    print("This means to separate most decay events, the detector system's time resolution must be significantly shorter than 1 ms.")

    # 2. Analyze the effect of the Distance between Detectors
    print("\nStep 2: Estimate the time-of-flight based on the distance.")
    # The electron from Bi-207 decay has ~1 MeV kinetic energy, so its speed is close to the speed of light.
    speed_of_light = 3.0e8  # meters per second
    distance_to_detector = 0.5 # meters (source is in the middle)
    time_of_flight_s = distance_to_detector / speed_of_light
    time_of_flight_ns = time_of_flight_s * 1e9 # convert to nanoseconds

    print(f"The time for an electron to travel 0.5 m to a detector is approximately {time_of_flight_ns:.2f} nanoseconds (ns).")
    
    # 3. Compare the time scales
    print("\nStep 3: Compare the time scales.")
    print(f"Time between decays (from Activity):  {avg_time_between_decays_s * 1e6:,.0f} microseconds (µs)")
    print(f"Time of flight (from Distance):       {time_of_flight_s * 1e6:.4f} microseconds (µs)")
    print("\nThe time between decays is millions of times longer than the electron's time-of-flight.")

    # 4. Conclusion
    print("\n--- Conclusion ---")
    print("The primary challenge in resolving individual decay events is determined by how frequently they occur.")
    print("The activity of the source (1000 decays/second) sets this fundamental rate.")
    print("The other factors, like time-of-flight (~1.7 ns), are orders of magnitude smaller and do not set the main requirement for separating one decay from the next.")
    print("\nTherefore, the dominant factor is the measured activity of the source.")


if __name__ == '__main__':
    analyze_time_resolution()
