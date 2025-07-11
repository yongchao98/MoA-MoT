import math

def calculate_molecular_speed():
    """
    This function demonstrates that molecular motion, the driving force of the
    experiment, is dependent on temperature.
    """
    # Constants
    R = 8.314  # Ideal gas constant in J/(mol·K)
    MOLAR_MASS_N2 = 0.028  # Molar mass of Nitrogen (N2) in kg/mol

    print("The experiment requires gas molecules to be in motion to pass through the one-way door.")
    print("This motion is a direct result of the gas's temperature.")
    print("We can calculate the typical speed of a gas molecule (root-mean-square speed) to demonstrate this.")
    print("-" * 50)

    # --- Case 1: A non-zero temperature (e.g., room temperature) ---
    temp_k = 298.15  # Room temperature in Kelvin
    
    # Calculation: v_rms = sqrt(3 * R * T / M)
    v_rms = math.sqrt((3 * R * temp_k) / MOLAR_MASS_N2)

    print(f"At a temperature of {temp_k} K:")
    print("The equation for rms speed is v_rms = sqrt((3 * R * T) / M)")
    print(f"v_rms = sqrt((3 * {R:.3f} * {temp_k}) / {MOLAR_MASS_N2})")
    print(f"The resulting molecular speed is approximately {v_rms:.2f} m/s.")
    print("This motion allows molecules to move between chambers.")
    print("-" * 50)

    # --- Case 2: Absolute zero ---
    temp_k_zero = 0.0
    
    # Calculation at absolute zero
    v_rms_zero = math.sqrt((3 * R * temp_k_zero) / MOLAR_MASS_N2)
    
    print(f"At a temperature of {temp_k_zero} K (absolute zero):")
    print(f"v_rms = sqrt((3 * {R:.3f} * {temp_k_zero}) / {MOLAR_MASS_N2})")
    print(f"The resulting molecular speed is {v_rms_zero:.2f} m/s.")
    print("At absolute zero, there is no molecular motion, so the experiment cannot proceed.")
    print("\nConclusion: Temperature is the required experimental parameter.")

calculate_molecular_speed()