import math

def analyze_time_resolution():
    """
    Calculates and compares the key time scales in the experiment to determine
    the dominant factor for the time resolution requirement.
    """

    # --- Given Parameters ---
    # The activity of the source in decays per second (1 kBq = 1000 Bq)
    activity_bq = 1000.0

    # The distance from the central source to one detector in meters
    distance_m = 0.5

    # --- Constants and Electron Properties ---
    # Speed of light in m/s
    C_M_S = 299792458
    # Approximate kinetic energy of a Bi-207 conversion electron in MeV
    KE_MEV = 1.0
    # Electron rest mass energy in MeV
    M0C2_MEV = 0.511

    # --- Calculation Step 1: Time scale from source activity ---
    # The average time between random decay events is the inverse of the activity rate.
    avg_time_between_decays_s = 1.0 / activity_bq

    # --- Calculation Step 2: Time scale from distance (Time of Flight) ---
    # Calculate the electron's speed. Since KE is comparable to rest mass, we use relativity.
    # Total Energy E = KE + m0c^2
    total_energy_mev = KE_MEV + M0C2_MEV
    # Relativistic factor gamma = E / m0c^2
    gamma = total_energy_mev / M0C2_MEV
    # Electron speed v = c * sqrt(1 - 1/gamma^2)
    electron_speed_m_s = C_M_S * math.sqrt(1 - (1 / gamma**2))
    # Time of flight = distance / speed
    time_of_flight_s = distance_m / electron_speed_m_s

    # --- Output and Conclusion ---
    print("Analysis of Time Resolution Requirements:")
    print("-" * 40)
    print("To individually measure electrons, the system's time resolution must be better")
    print("than the average time separating two distinct decay events.")
    print("\nComparing the relevant time scales:")
    
    print("\n1. Time Between Decays (from Activity):")
    print(f"   Source Activity = {int(activity_bq)} decays/second")
    print(f"   Average Time Between Decays = 1 / {int(activity_bq)} = {avg_time_between_decays_s:.3f} seconds")
    print(f"   This is equivalent to {avg_time_between_decays_s * 1000:.1f} milliseconds (ms).")

    print("\n2. Time of Flight (from Distance):")
    print(f"   Distance to Detector = {distance_m} meters")
    print(f"   Electron Time of Flight = {time_of_flight_s * 1e9:.2f} nanoseconds (ns).")

    print("\nConclusion:")
    print(f"The average time between decays is ~{avg_time_between_decays_s * 1000:.1f} ms, while the time of flight is ~{time_of_flight_s * 1e9:.2f} ns.")
    print("A millisecond is one million times longer than a nanosecond.")
    print("The time between decay events is vastly longer than the electron's travel time.")
    print("Therefore, the dominant factor that limits our ability to distinguish one event")
    print("from the next is the rate at which they occur, which is the source activity.")

analyze_time_resolution()