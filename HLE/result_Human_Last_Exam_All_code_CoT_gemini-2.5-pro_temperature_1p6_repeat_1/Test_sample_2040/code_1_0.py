import math

def calculate_spdc_wavelengths(pump_wavelength_nm):
    """
    Calculates and prints wavelengths for hypothetical SPDC scenarios.

    Args:
        pump_wavelength_nm (float): The wavelength of the pump photon in nanometers.
    """
    print(f"--- Assuming SPDC in a hypothetical material ---")
    print(f"Pump Wavelength (λ_pump): {pump_wavelength_nm} nm\n")
    print("The governing principle is energy conservation, which for wavelengths is:")
    print("1/λ_pump = 1/λ_signal + 1/λ_idler\n")

    # --- Case 1: Degenerate SPDC ---
    # The pump photon splits into two identical photons.
    # 1/λ_pump = 1/λ_signal + 1/λ_signal = 2/λ_signal
    # So, λ_signal = 2 * λ_pump
    degenerate_output_wavelength_nm = 2 * pump_wavelength_nm
    print("--- Case 1: Degenerate SPDC (photons are identical) ---")
    print(f"Calculation: λ_signal = λ_idler = 2 * λ_pump")
    # Print the equation with numbers
    print(f"Result: {degenerate_output_wavelength_nm:.2f} nm = 2 * {pump_wavelength_nm} nm")
    # Verify the result with the original formula
    # Note: 1/pump = 1/(2*pump) + 1/(2*pump) = 2/(2*pump) = 1/pump
    print(f"Final Equation: 1 / {pump_wavelength_nm} = 1 / {degenerate_output_wavelength_nm:.2f} + 1 / {degenerate_output_wavelength_nm:.2f}\n")


    # --- Case 2: Non-Degenerate SPDC ---
    # The pump photon splits into two different photons. We must define one to find the other.
    # Let's assume we detect a signal photon at a specific wavelength.
    # 1/λ_idler = 1/λ_pump - 1/λ_signal
    # λ_idler = 1 / (1/λ_pump - 1/λ_signal)
    signal_wavelength_nm = 810.0  # A common Ti:Sapphire laser wavelength
    
    # Check if this is physically possible (signal wavelength must be > pump wavelength)
    if signal_wavelength_nm <= pump_wavelength_nm:
        print(f"--- Case 2: Non-Degenerate SPDC (Error) ---")
        print(f"Signal wavelength ({signal_wavelength_nm} nm) must be longer than pump wavelength ({pump_wavelength_nm} nm).")
        return
        
    idler_inv_wavelength = (1 / pump_wavelength_nm) - (1 / signal_wavelength_nm)
    idler_wavelength_nm = 1 / idler_inv_wavelength

    print("--- Case 2: Non-Degenerate SPDC (photons are different) ---")
    print(f"Calculation: Given a signal photon of {signal_wavelength_nm} nm, we calculate the idler photon.")
    print("Formula: 1/λ_idler = 1/λ_pump - 1/λ_signal")
    print(f"Resulting Idler Wavelength (λ_idler): {idler_wavelength_nm:.2f} nm")
    # Print the final equation with all numbers
    print(f"Final Equation: 1 / {pump_wavelength_nm} = 1 / {signal_wavelength_nm:.2f} + 1 / {idler_wavelength_nm:.2f}")


# Let's use a common green laser pointer wavelength for the pump
pump_laser_wavelength = 532.0
calculate_spdc_wavelengths(pump_laser_wavelength)