import math

def analyze_timescales():
    """
    Calculates and compares the different timescales in the Bi-207 decay experiment
    to determine the dominant factor for time resolution.
    """
    # Given parameters
    activity_bq = 1000.0  # Source activity in Bq (1 kBq)
    distance_m = 0.5  # Distance from the source to one detector
    electron_ke_kev = 976.0  # A typical kinetic energy for a Bi-207 conversion electron in keV

    # Physical constants
    electron_rest_mass_kev = 511.0  # Electron rest mass in keV/c^2
    speed_of_light_m_s = 299792458.0  # Speed of light in m/s
    correlated_emission_halflife_s = 129e-12  # Half-life of Pb-207m state is ~129 ps

    # --- Step 1: Calculate the average time between decays from activity ---
    avg_time_between_decays_s = 1.0 / activity_bq

    # --- Step 2: Calculate the electron time-of-flight from distance and energy ---
    # Relativistic calculation for electron velocity
    gamma = (electron_ke_kev / electron_rest_mass_kev) + 1.0
    beta = math.sqrt(1.0 - 1.0 / (gamma**2))  # beta = v/c
    electron_velocity_m_s = beta * speed_of_light_m_s
    time_of_flight_s = distance_m / electron_velocity_m_s

    # --- Step 3: Present the comparison ---
    print("Analysis of Characteristic Timescales:")
    print("-" * 40)

    # Timescale from Source Activity
    print(f"1. Timescale from Source Activity (1000 Bq):")
    print(f"   The average time between individual decays is 1 / {int(activity_bq)} Hz = {avg_time_between_decays_s * 1000:.1f} milliseconds.")

    # Timescale from Distance (Time-of-Flight)
    print(f"\n2. Timescale from Detector Distance (0.5 m):")
    print(f"   The time-of-flight for a {int(electron_ke_kev)} keV electron is {time_of_flight_s * 1e9:.2f} nanoseconds.")

    # Timescale from Correlated Emissions
    print(f"\n3. Timescale from Correlated Emissions:")
    print(f"   The time between cascade emissions is on the order of {correlated_emission_halflife_s * 1e12:.0f} picoseconds.")
    
    print("-" * 40)
    print("\nConclusion:")
    print("To 'individually measure the electrons', the system must distinguish electrons from different decays.")
    print(f"The average time between these decays ({avg_time_between_decays_s * 1000:.1f} ms) is much longer than any other timescale in the system.")
    print("Therefore, the system's time resolution must be fast enough to handle events arriving, on average, every millisecond.")
    print("This requirement is set by the rate of events, which is the source activity.")

if __name__ == '__main__':
    analyze_timescales()