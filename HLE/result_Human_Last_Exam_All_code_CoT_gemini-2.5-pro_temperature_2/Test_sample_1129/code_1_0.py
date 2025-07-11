import numpy as np

def calculate_fusion_atoms():
    """
    Calculates the minimum number of Ti-50 atoms for a hypothetical
    sustained fusion reaction based on the Lawson criterion.
    """
    # Step 1: Define constants and given values
    
    # Boltzmann constant in eV/K
    BOLTZMANN_CONSTANT_EV_K = 8.617333e-5
    
    # Lawson Criterion value (L) for D-T fusion, used as a proxy.
    # Units: keV * s / m^3
    LAWSON_CRITERION_DT = 3e21
    
    # Temperature in Celsius (given)
    temp_celsius = 20.0
    
    # Confinement time in seconds (given)
    confinement_time_s = 1.0
    
    # Side length of the reaction chamber in cm (given)
    side_length_cm = 10.0
    
    # --- Calculations ---
    
    # Step 2: Convert temperature to Kelvin and then to keV
    temp_kelvin = temp_celsius + 273.15
    # Energy in eV = k_B * T
    temp_ev = BOLTZMANN_CONSTANT_EV_K * temp_kelvin
    # Temperature (T) in keV
    temp_kev = temp_ev / 1000.0
    
    # Step 3: Calculate required particle density (n) using Lawson criterion
    # n * T * τ ≥ L  =>  n = L / (T * τ)
    # Units: (keV*s/m^3) / (keV * s) = 1/m^3
    particle_density_m3 = LAWSON_CRITERION_DT / (temp_kev * confinement_time_s)
    
    # Step 4: Calculate reaction chamber volume (V)
    # Convert cm to m
    side_length_m = side_length_cm / 100.0
    # Volume in m^3
    volume_m3 = side_length_m ** 3
    
    # Step 5: Calculate total number of atoms (N)
    # N = n * V
    total_atoms = particle_density_m3 * volume_m3
    
    # --- Output the results ---
    
    print("This calculation assumes a hypothetical room-temperature Ti-50 fusion reaction.")
    print("The Lawson Criterion value for D-T fusion is used as a proxy.\n")
    print("--- Equation and Values ---")
    print(f"Number of Atoms = (Lawson Criterion / (Temperature * Confinement Time)) * Volume")
    
    # I will format the numbers to show them clearly in the equation
    lawson_str = f"{LAWSON_CRITERION_DT:.1e}"
    temp_str = f"{temp_kev:.4e}"
    time_str = f"{confinement_time_s}"
    vol_str = f"{volume_m3}"
    
    print(f"Number of Atoms = ({lawson_str} / ({temp_str} * {time_str})) * {vol_str}")
    print("\n--- Result ---")
    print(f"The minimum number of titanium-50 atoms required is: {total_atoms:.4e}")

# Run the calculation
calculate_fusion_atoms()