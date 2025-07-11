import math

def calculate_timescales():
    """
    Calculates and compares the time resolution requirements set by
    the source activity and the detector distance.
    """
    # --- Constants ---
    C_M_PER_S = 299792458  # Speed of light in m/s
    M_E_C2_MEV = 0.511      # Electron rest mass energy in MeV

    # --- Experimental Parameters ---
    ACTIVITY_Bq = 1000  # 1 kBq = 1000 decays per second
    DISTANCE_M = 1.0      # Distance between detectors
    SOURCE_TO_DETECTOR_M = DISTANCE_M / 2

    # --- 1. Timescale from Source Activity ---
    # The activity determines the average time between decay events.
    # To avoid accidental coincidences, the time resolution must be much smaller than this.
    avg_time_between_decays_s = 1.0 / ACTIVITY_Bq
    
    print("--- Analysis of Timescales ---")
    print(f"Source Activity: {ACTIVITY_Bq} decays/second")
    print(f"Average time between decays: {avg_time_between_decays_s * 1e3:.2f} ms ({avg_time_between_decays_s * 1e9:.2f} ns)")
    print("This means the requirement to avoid accidental coincidences is that the time window should be MUCH LESS than 1 ms.\n")

    # --- 2. Timescale from Distance (Time-of-Flight) ---
    # The distance determines the time-of-flight (ToF) of the electrons.
    # In a coincidence experiment, we need to resolve the ToF difference
    # between two correlated electrons of different energies.
    # We'll use two prominent conversion electrons from Bi-207 decay as an example.
    e1_ke_mev = 0.976  # Energy of electron 1 (976 keV)
    e2_ke_mev = 0.482  # Energy of electron 2 (482 keV)

    def get_tof_ns(ke_mev):
        """Calculates the time of flight in nanoseconds for an electron."""
        # Relativistic gamma factor: KE = (gamma - 1) * m_e * c^2
        gamma = (ke_mev / M_E_C2_MEV) + 1
        # Velocity: v = c * sqrt(1 - 1/gamma^2)
        beta_sq = 1 - (1 / (gamma**2))
        velocity_m_per_s = C_M_PER_S * math.sqrt(beta_sq)
        # Time of flight: t = d / v
        tof_s = SOURCE_TO_DETECTOR_M / velocity_m_per_s
        return tof_s * 1e9 # convert to nanoseconds

    tof1_ns = get_tof_ns(e1_ke_mev)
    tof2_ns = get_tof_ns(e2_ke_mev)
    tof_difference_ns = abs(tof1_ns - tof2_ns)

    print("--- Time-of-Flight Calculation ---")
    print(f"Distance from source to each detector: {SOURCE_TO_DETECTOR_M} m")
    print(f"Time of flight for {e1_ke_mev*1000:.0f} keV electron: {tof1_ns:.2f} ns")
    print(f"Time of flight for {e2_ke_mev*1000:.0f} keV electron: {tof2_ns:.2f} ns")
    print(f"Time difference of arrival at detectors: {tof_difference_ns:.2f} ns")
    print("This means the time resolution must be on the order of nanoseconds (or better) to distinguish the arrival of these two electrons.\n")

    # --- 3. Conclusion ---
    print("--- Conclusion ---")
    print(f"The requirement from activity is to have a resolution much less than {avg_time_between_decays_s * 1e3:.0f} ms.")
    print(f"The requirement from the detector distance is to have a resolution around {tof_difference_ns:.2f} ns.")
    print(f"Since {tof_difference_ns:.2f} ns is a much shorter (more stringent) time requirement than {avg_time_between_decays_s * 1e6:.0f} ns (1 ms), the dominant factor setting the minimum time resolution is the distance between the detectors.")

calculate_timescales()