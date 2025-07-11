import math

def solve_timescales():
    """
    Calculates and compares the timescales set by source activity and particle time-of-flight.
    """
    # --- GIVEN PARAMETERS ---
    # Source Activity
    activity_bq = 1000  # 1 kBq = 1000 decays/second

    # Experiment Geometry and Particle Physics
    distance_m = 0.5    # Distance from source to a single detector
    ke_ev = 1.0e6       # Kinetic energy of electron in eV (1 MeV)

    # --- PHYSICAL CONSTANTS ---
    c_ms = 299792458                  # Speed of light in m/s
    m_e_ev = 510998.9                 # Electron rest mass in eV/c^2

    # --- CALCULATIONS ---

    # 1. Timescale from Source Activity
    # The average time between random decay events is the inverse of the decay rate (activity).
    avg_time_between_decays_s = 1.0 / activity_bq

    # 2. Timescale from Distance (Time of Flight)
    # To find the electron's speed, we use relativistic kinematics.
    # Total Energy E = Kinetic Energy (KE) + Rest Mass Energy (m_e*c^2)
    total_energy_ev = ke_ev + m_e_ev
    # The Lorentz factor, gamma = Total Energy / Rest Mass Energy
    gamma = total_energy_ev / m_e_ev
    # The speed v is related to gamma by: gamma = 1 / sqrt(1 - v^2/c^2)
    # Rearranging for v gives: v = c * sqrt(1 - 1/gamma^2)
    electron_velocity_ms = c_ms * math.sqrt(1 - (1 / gamma**2))
    # Time of Flight = Distance / Velocity
    time_of_flight_s = distance_m / electron_velocity_ms

    # --- OUTPUT ---
    print("To determine the dominant factor for time resolution, we compare the key timescales.")
    print("-" * 70)
    
    print(f"Factor 1: Timescale set by the Source Activity ({activity_bq} Bq)")
    print(f"The average time between two separate decay events is:")
    print(f"  T_decay = 1 / {activity_bq} Hz = {avg_time_between_decays_s} s")
    print(f"  Which is {avg_time_between_decays_s * 1e3:.1f} milliseconds or {avg_time_between_decays_s * 1e9:,.0f} nanoseconds.")
    
    print("\n" + "-" * 70)
    
    print(f"Factor 2: Timescale set by the Detector Distance ({distance_m * 2} m apart)")
    print(f"The time for a 1 MeV electron to travel {distance_m} m to a detector is:")
    print(f"  t_flight = {distance_m} m / {electron_velocity_ms:.2e} m/s = {time_of_flight_s:.2e} s")
    print(f"  Which is {time_of_flight_s * 1e9:.2f} nanoseconds.")

    print("\n" + "-" * 70)
    print("Conclusion:")
    print(f"The average time between decays ({avg_time_between_decays_s * 1e6:.1f} microseconds) is much larger than the electron's flight time ({time_of_flight_s * 1e9:.2f} nanoseconds).")
    print("The system must be able to distinguish events that are randomly spaced in time, with an average separation set by the activity.")
    print("Therefore, the dominant factor that sets the minimum time resolution requirement to separate one decay event from the next is the activity of the source.")

solve_timescales()
<<<D>>>