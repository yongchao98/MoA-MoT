import math

def calculate_rms_speed():
    """
    Calculates the root-mean-square (RMS) speed of gas molecules.
    This demonstrates that molecular speed is dependent on temperature, which
    is the required parameter for the Maxwell's demon apparatus to function.
    """
    # --- Constants ---
    # T: Temperature in Kelvin (e.g., room temperature 25°C = 298.15 K)
    T = 298.15
    # k_B: Boltzmann constant in J/K
    k_B = 1.380649e-23
    # M_N2: Molar mass of Nitrogen (N2) in kg/mol
    M_N2 = 0.0280134
    # N_A: Avogadro's number in mol^-1
    N_A = 6.02214076e23

    # --- Calculation ---
    # m: Mass of a single nitrogen molecule in kg
    m = M_N2 / N_A
    
    # v_rms: Root-mean-square speed in m/s
    v_rms = math.sqrt((3 * k_B * T) / m)

    # --- Output ---
    print("This script calculates the speed of gas molecules based on temperature.")
    print("In the thought experiment, this speed is necessary for molecules to cross the one-way door.")
    print(f"If Temperature (T) were 0, the speed would be 0, and no gas would be trapped.\n")
    print(f"Calculation for Nitrogen gas at {T} K (room temperature):")
    print("-" * 50)
    print(f"Formula: v_rms = sqrt(3 * k_B * T / m)\n")
    print("Plugging in the values:")
    # The final equation with all numbers is printed below
    print(f"v_rms = sqrt(3 * {k_B:.6e} J/K * {T} K / {m:.6e} kg)")
    print(f"\nCalculated RMS Speed of a Nitrogen Molecule: {v_rms:.2f} m/s")

calculate_rms_speed()