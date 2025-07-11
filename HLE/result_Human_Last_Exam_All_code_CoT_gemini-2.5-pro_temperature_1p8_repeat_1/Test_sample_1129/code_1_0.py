import scipy.constants

def calculate_fusion_atoms():
    """
    Calculates the minimum number of Titanium-50 atoms for a hypothetical
    sustained fusion reaction based on the Lawson criterion under the specified
    room temperature conditions.
    """

    # --- Step 1: Define constants and given parameters ---
    # Lawson criterion for D-T ignition (used for this hypothetical reaction)
    # Units: keV * s * m^-3
    lawson_constant_C = 5e21

    # Confinement time (tau), in seconds
    confinement_time_tau = 1.0  # s

    # Room temperature, in Celsius
    room_temp_celsius = 20.0  # °C

    # Reaction chamber side length, in cm
    chamber_side_length_cm = 10.0  # cm

    # --- Step 2: Convert units for calculation ---

    # Convert temperature from Celsius to Kelvin
    temp_k = room_temp_celsius + scipy.constants.zero_Celsius

    # Convert temperature from Kelvin to keV
    # 1 eV / k_B = 11604.5 K. So 1 keV = 11604500 K.
    K_PER_KEV = 1.16045e7
    temp_kev = temp_k / K_PER_KEV

    # Convert chamber side length from cm to meters
    chamber_side_length_m = chamber_side_length_cm / 100.0

    print("--- Problem Setup ---")
    print(f"Lawson Criterion Constant (C): {lawson_constant_C:.1e} keV*s/m^3")
    print(f"Confinement Time (τ): {confinement_time_tau} s")
    print(f"Temperature (T): {room_temp_celsius}°C = {temp_k:.2f} K = {temp_kev:.5e} keV")
    print(f"Chamber Side Length: {chamber_side_length_cm} cm = {chamber_side_length_m} m\n")

    # --- Step 3: Calculate the required particle density (n) ---
    # From Lawson Criterion: n * τ * T >= C  =>  n >= C / (τ * T)
    required_density_n = lawson_constant_C / (confinement_time_tau * temp_kev)

    print("--- Calculation Steps ---")
    print(f"1. Required Particle Density (n) = C / (τ * T)")
    print(f"   n = {lawson_constant_C:.1e} / ({confinement_time_tau} * {temp_kev:.5e})")
    print(f"   n ≈ {required_density_n:.4e} particles/m^3\n")


    # --- Step 4: Calculate the volume of the chamber (V) ---
    chamber_volume_m3 = chamber_side_length_m ** 3

    print(f"2. Chamber Volume (V) = (Side Length)^3")
    print(f"   V = ({chamber_side_length_m})^3 = {chamber_volume_m3} m^3\n")


    # --- Step 5: Calculate the total number of atoms (N) ---
    # N = n * V
    total_atoms_N = required_density_n * chamber_volume_m3

    print(f"3. Total Number of Atoms (N) = n * V")
    print("   Final Equation:")
    print(f"   N = ({lawson_constant_C} / ({confinement_time_tau} * {temp_kev})) * {chamber_volume_m3}")
    print("\n--- Final Result ---")
    print(f"The minimum number of Titanium-50 atoms required is: {total_atoms_N:.4e}")

if __name__ == '__main__':
    calculate_fusion_atoms()
    # The final answer from the calculation is approximately 1.9793e+23
    final_answer = 1.979280043836791e+23
    print(f"\n<<<1.9793e+23>>>")
