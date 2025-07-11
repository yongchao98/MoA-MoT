def calculate_spdc_idler_wavelength(pump_wavelength_nm, signal_wavelength_nm):
    """
    In SPDC, energy is conserved: E_pump = E_signal + E_idler.
    Since energy E is inversely proportional to wavelength λ (E = hc/λ),
    the conservation law can be written as:
    1/λ_pump = 1/λ_signal + 1/λ_idler

    This function calculates the idler wavelength based on this principle.

    Args:
        pump_wavelength_nm (float): The wavelength of the pump photon in nm.
        signal_wavelength_nm (float): The wavelength of the signal photon in nm.

    Returns:
        float: The wavelength of the idler photon in nm, or None if inputs are invalid.
    """
    if pump_wavelength_nm <= 0 or signal_wavelength_nm <= 0:
        print("Error: Wavelengths must be positive values.")
        return None
    
    # The signal photon must have lower energy (longer wavelength) than the pump photon.
    if signal_wavelength_nm <= pump_wavelength_nm:
        print(f"Error: Signal wavelength ({signal_wavelength_nm} nm) must be longer than pump wavelength ({pump_wavelength_nm} nm).")
        return None

    # Calculate 1/λ_idler from 1/λ_idler = 1/λ_pump - 1/λ_signal
    try:
        inv_idler_wavelength = (1 / pump_wavelength_nm) - (1 / signal_wavelength_nm)
        if inv_idler_wavelength <= 1e-9: # Effectively zero or negative
             return None
        idler_wavelength_nm = 1 / inv_idler_wavelength
        return idler_wavelength_nm
    except ZeroDivisionError:
        return None

# --- Example Calculation ---
# This is a hypothetical calculation to demonstrate the physical principle.
# We assume a common green laser for the pump and a standard near-infrared
# wavelength for the signal photon.

# Input values
pump_lambda = 532.0  # Wavelength of a common green laser pointer (nm)
signal_lambda = 810.0 # Wavelength of a common Ti:Sapphire laser (nm)

# Perform the calculation
idler_lambda = calculate_spdc_idler_wavelength(pump_lambda, signal_lambda)

# Print the results in the requested format
if idler_lambda is not None:
    print("This script demonstrates the energy conservation principle in a hypothetical SPDC process.")
    print(f"Given a pump wavelength of {pump_lambda:.1f} nm and a signal wavelength of {signal_lambda:.1f} nm, the idler wavelength is calculated.")
    print("\n--- Result ---")
    print(f"The final equation showing the down-conversion is:")
    print(f"{pump_lambda:.1f} nm (pump) -> {signal_lambda:.1f} nm (signal) + {idler_lambda:.1f} nm (idler)")
else:
    print("Could not calculate a valid idler wavelength with the provided inputs.")
