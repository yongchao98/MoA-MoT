# -*- coding: utf-8 -*-

def calculate_spdc_wavelengths():
    """
    This script demonstrates the principle of energy conservation in
    Spontaneous Parametric Down-Conversion (SPDC), a process not expected
    in free-standing boron nanosheets due to their crystal symmetry.
    """
    # Define pump and signal wavelengths in nanometers (nm).
    # A 405 nm laser is a common pump source.
    # We'll use the degenerate case where the pump photon splits into two
    # identical photons, so the signal/idler wavelength is double the pump.
    pump_wavelength = 405.0
    signal_wavelength = 810.0

    # According to the law of energy conservation for SPDC:
    # E_pump = E_signal + E_idler
    # Since E = hc/λ, this becomes: hc/λ_p = hc/λ_s + hc/λ_i
    # which simplifies to: 1/λ_p = 1/λ_s + 1/λ_i

    # We can rearrange to solve for the idler wavelength:
    # 1/λ_i = 1/λ_p - 1/λ_s
    try:
        # Avoid division by zero and ensure the physical condition that
        # signal wavelength is greater than pump wavelength.
        if pump_wavelength <= 0 or signal_wavelength <= pump_wavelength:
            raise ValueError("Wavelengths must be positive, and signal wavelength must be longer than pump wavelength.")
            
        inv_idler_wavelength = (1 / pump_wavelength) - (1 / signal_wavelength)
        idler_wavelength = 1 / inv_idler_wavelength
    
        print("This script demonstrates energy conservation in Spontaneous Parametric Down-Conversion (SPDC).")
        print("\nBased on the energy conservation equation:")
        print("1 / (Pump Wavelength) = 1 / (Signal Wavelength) + 1 / (Idler Wavelength)\n")
        print(f"For a pump photon at {pump_wavelength:.1f} nm and a detected signal photon at {signal_wavelength:.1f} nm,")
        print(f"the corresponding idler photon must have a wavelength of {idler_wavelength:.1f} nm.\n")

        print("The full equation with these values is:")
        print(f"1 / {pump_wavelength:.1f} = 1 / {signal_wavelength:.1f} + 1 / {idler_wavelength:.1f}")

    except (ValueError, ZeroDivisionError) as e:
        print(f"An error occurred during calculation: {e}")

if __name__ == '__main__':
    calculate_spdc_wavelengths()
