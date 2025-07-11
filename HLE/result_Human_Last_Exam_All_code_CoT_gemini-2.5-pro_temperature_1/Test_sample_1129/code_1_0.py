import math

def calculate_fusion_atoms():
    """
    Calculates the number of atoms in a given volume at a specific
    temperature and pressure using the Ideal Gas Law.
    """
    # Given parameters
    side_length_cm = 10.0
    temperature_c = 20.0
    pressure_pa = 101325.0

    # Constants
    k_boltzmann = 1.380649e-23  # Boltzmann constant in J/K

    # --- Step 1: Convert all units to SI ---
    # Convert side length from cm to meters
    side_length_m = side_length_cm / 100.0
    # Calculate volume in cubic meters
    volume_m3 = side_length_m ** 3
    # Convert temperature from Celsius to Kelvin
    temperature_k = temperature_c + 273.15

    # --- Step 2: Apply the Ideal Gas Law to find N ---
    # The Ideal Gas Law is PV = N * k_B * T
    # We solve for N: N = (P * V) / (k_B * T)
    num_atoms = (pressure_pa * volume_m3) / (k_boltzmann * temperature_k)

    # --- Step 3: Print the results ---
    print("To find the minimum number of atoms required for the reaction under the specified conditions, we use the Ideal Gas Law.")
    print("\nGiven conditions:")
    print(f"Pressure (P): {pressure_pa} Pa")
    print(f"Volume (V): {volume_m3:.3f} m^3 (from a {side_length_cm} cm cube)")
    print(f"Temperature (T): {temperature_k:.2f} K (from {temperature_c}°C)")
    print(f"Boltzmann Constant (k_B): {k_boltzmann:.6e} J/K")

    print("\nIdeal Gas Law formula to find the number of atoms (N):")
    print("N = (P * V) / (k_B * T)")

    print("\nSubstituting the values into the equation:")
    # The final equation with each number explicitly shown
    print(f"N = ({pressure_pa} * {volume_m3:.3f}) / ({k_boltzmann:.6e} * {temperature_k:.2f})")
    
    numerator = pressure_pa * volume_m3
    denominator = k_boltzmann * temperature_k
    
    print(f"N = {numerator:.3f} / {denominator:.6e}")
    print(f"N = {num_atoms:.4e}")

    print(f"\nTherefore, the minimum number of titanium-50 atoms required is approximately {num_atoms:.4e}.")
    
    # Final answer in the specified format
    print(f"\n<<<{num_atoms}>>>")

if __name__ == '__main__':
    calculate_fusion_atoms()