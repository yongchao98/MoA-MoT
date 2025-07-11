import sys

def calculate_spdc_energy_conservation(pump_wavelength_nm, signal_wavelength_nm=None):
    """
    Demonstrates energy conservation in Spontaneous Parametric Down-Conversion (SPDC).

    The energy conservation law is E_pump = E_signal + E_idler.
    Since E = hc/λ, this is equivalent to: 1/λ_pump = 1/λ_signal + 1/λ_idler.

    This function showcases two cases:
    1. Degenerate SPDC: The signal and idler photons have the same wavelength.
    2. Non-degenerate SPDC: The signal and idler photons have different wavelengths.

    Args:
        pump_wavelength_nm (float): Wavelength of the pump photon in nanometers.
        signal_wavelength_nm (float, optional): Wavelength of the signal photon in nm for the non-degenerate case.
                                                 If not provided, only the degenerate case is shown.
    """
    print(f"Analyzing SPDC for a pump laser with wavelength: {pump_wavelength_nm} nm\n")

    # --- Case 1: Degenerate SPDC ---
    # The pump photon splits into two identical daughter photons.
    degenerate_photon_wavelength = 2 * pump_wavelength_nm
    print("--- 1. Degenerate SPDC Example ---")
    print("The pump photon splits into two photons of equal energy and wavelength.")
    print("The energy conservation equation is: 1/λ_pump = 1/λ_signal + 1/λ_idler")
    print(f"The equation with values is: 1 / {pump_wavelength_nm} = 1 / {degenerate_photon_wavelength:.2f} + 1 / {degenerate_photon_wavelength:.2f}")
    print("-" * 35)

    # --- Case 2: Non-Degenerate SPDC ---
    if signal_wavelength_nm is not None:
        print("\n--- 2. Non-Degenerate SPDC Example ---")
        # The signal wavelength must be greater than the pump wavelength for energy to be conserved.
        if signal_wavelength_nm <= pump_wavelength_nm:
            print(f"Error: Signal wavelength ({signal_wavelength_nm} nm) must be longer than pump wavelength ({pump_wavelength_nm} nm).", file=sys.stderr)
            return

        # Calculate idler wavelength from: 1/λ_idler = 1/λ_pump - 1/λ_signal
        inv_idler_wavelength = (1 / pump_wavelength_nm) - (1 / signal_wavelength_nm)

        if inv_idler_wavelength <= 0:
             print(f"Error: Calculation resulted in a non-physical idler wavelength.", file=sys.stderr)
             return

        idler_wavelength_nm = 1 / inv_idler_wavelength
        print(f"Given a detected signal photon of {signal_wavelength_nm} nm, the idler photon is calculated.")
        print("The energy conservation equation is: 1/λ_pump = 1/λ_signal + 1/λ_idler")
        print(f"The equation with values is: 1 / {pump_wavelength_nm} = 1 / {signal_wavelength_nm} + 1 / {idler_wavelength_nm:.2f}")
        print("-" * 35)


if __name__ == '__main__':
    # --- Execute the calculations with example values ---
    # Common pump laser (e.g., from a Blu-ray diode)
    PUMP_LAMBDA = 405.0 # nm
    # A possible signal wavelength for the non-degenerate case
    SIGNAL_LAMBDA = 780.0 # nm

    calculate_spdc_energy_conservation(PUMP_LAMBDA, SIGNAL_LAMBDA)
