import math

def find_dominant_factor():
    """
    Calculates and compares the timescales related to source activity and detector
    distance to find the dominant factor for time resolution.
    """
    # --- Input Parameters from the problem statement ---
    activity_Bq = 1000.0  # 1 kBq = 1000 decays/second
    distance_to_detector_m = 0.5 # Source is in the middle of 1m span
    electron_KE_MeV = 1.0  # A typical energy for Bi-207 conversion electrons

    # --- Physical Constants ---
    c_m_per_s = 299792458.0      # Speed of light in m/s
    electron_rest_mass_MeV = 0.511 # Electron rest mass energy in MeV

    # --- Step 1: Calculate the average time between decays from the source activity ---
    print("--- Analysis of Factor (D): Source Activity ---")
    avg_time_between_decays_s = 1.0 / activity_Bq
    print("The measured activity determines the average rate of electron emissions.")
    print("The average time between consecutive decays is the inverse of the activity.")
    print(f"Calculation: 1 / {int(activity_Bq)} decays/sec = {avg_time_between_decays_s:.3f} seconds")
    print(f"(This is equivalent to {avg_time_between_decays_s * 1000:.0f} milliseconds)\n")

    # --- Step 2: Calculate the electron's time-of-flight based on detector distance ---
    print("--- Analysis of Factor (A): Detector Distance ---")
    # Using relativistic kinetics to find the electron's velocity
    total_energy_MeV = electron_KE_MeV + electron_rest_mass_MeV
    gamma = total_energy_MeV / electron_rest_mass_MeV
    # beta = v/c = sqrt(1 - 1/gamma^2)
    beta = math.sqrt(1 - 1 / (gamma**2))
    electron_velocity_m_per_s = beta * c_m_per_s
    time_of_flight_s = distance_to_detector_m / electron_velocity_m_per_s
    
    print("The detector distance determines the electron's time-of-flight from source to detector.")
    print("First, we calculate the velocity of a 1.0 MeV electron.")
    print(f"Velocity Calculation: v = c * sqrt(1 - (m_e*c^2 / E_total)^2) = {electron_velocity_m_per_s:.3e} m/s")
    print("Next, we calculate the time-of-flight.")
    print(f"Calculation: {distance_to_detector_m} m / {electron_velocity_m_per_s:.3e} m/s = {time_of_flight_s:.3e} seconds")
    print(f"(This is equivalent to {time_of_flight_s * 1e9:.2f} nanoseconds)\n")
    
    # --- Step 3: Compare the timescales and conclude ---
    print("--- Conclusion ---")
    print(f"Average time between decays: {avg_time_between_decays_s * 1e3:.1f} ms")
    print(f"Electron time of flight:     {time_of_flight_s * 1e9:.2f} ns")
    print("\nThe system must be able to distinguish events that are, on average, 1 millisecond apart.")
    print("The time it takes for an electron to reach the detector (~1.8 nanoseconds) is about 500,000 times shorter than the average time between the events we want to measure individually.")
    print("Therefore, the dominant factor that sets the time resolution requirement is the rate at which events occur, which is determined by the source activity.")

find_dominant_factor()
<<<D>>>