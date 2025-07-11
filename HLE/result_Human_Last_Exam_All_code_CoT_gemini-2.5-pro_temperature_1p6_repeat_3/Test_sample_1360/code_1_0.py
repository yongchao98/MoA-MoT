import math

def analyze_time_resolution_factors():
    """
    Analyzes the factors affecting time resolution for detecting Bi-207 decay electrons.
    """

    # 1. Analyze the effect of source activity
    activity_kBq = 1
    activity_Bq = activity_kBq * 1000  # Convert kBq to Bq (decays per second)

    print("Step 1: Calculate the average time between decay events from the source activity.")
    # The average time interval between decays is the inverse of the activity rate.
    avg_time_between_decays_s = 1 / activity_Bq
    avg_time_between_decays_ms = avg_time_between_decays_s * 1000

    print(f"The source activity is {activity_Bq} decays/second.")
    print("The average time interval (Δt) between decays is calculated as: 1 / Activity")
    print(f"Δt = 1 / {activity_Bq} decays/second = {avg_time_between_decays_s:.3f} seconds")
    print(f"This is equal to {avg_time_between_decays_ms:.1f} milliseconds.\n")
    print("To distinguish individual decay events, the system's time resolution must be significantly better than this value.")

    # 2. Analyze the effect of the distance between detectors
    distance_m = 0.5  # Distance from source to one detector
    # Bi-207 emits conversion electrons with energies up to ~1 MeV.
    # An electron with 1 MeV kinetic energy is relativistic and travels at about 94% of the speed of light.
    speed_of_light_m_s = 3.0e8
    electron_velocity_m_s = 0.94 * speed_of_light_m_s

    print("Step 2: Calculate the electron's time-of-flight from the source to a detector.")
    # The time of flight (TOF) is distance / velocity.
    tof_s = distance_m / electron_velocity_m_s
    tof_ns = tof_s * 1e9

    print(f"The distance from the source to a detector is {distance_m} m.")
    print(f"The approximate speed of a high-energy electron is {electron_velocity_m_s:.2e} m/s.")
    print("The time-of-flight (TOF) is calculated as: Distance / Velocity")
    print(f"TOF = {distance_m} m / {electron_velocity_m_s:.2e} m/s = {tof_s:.2e} seconds")
    print(f"This is equal to {tof_ns:.2f} nanoseconds.\n")

    # 3. Compare the time scales
    print("Step 3: Compare the time scales.")
    print(f"Time between decays: {avg_time_between_decays_ms:.1f} ms")
    print(f"Electron time-of-flight: {tof_ns:.2f} ns")
    print(f"\nConclusion: The average time between decays ({avg_time_between_decays_ms:.1f} ms) is much larger than the electron's time-of-flight ({tof_ns:.2f} ns).")
    print("Therefore, the dominant factor that determines the necessary time resolution to measure individual decay events is the rate of those events, which is given by the source activity.")

if __name__ == '__main__':
    analyze_time_resolution_factors()