def calculate_spdc_wavelengths():
    """
    This function illustrates the principle of energy conservation in
    Spontaneous Parametric Down-Conversion (SPDC).

    SPDC is a process where one high-energy 'pump' photon splits into two
    lower-energy photons: a 'signal' and an 'idler'. Energy conservation
    dictates that the energy of the pump equals the sum of the energies
    of the signal and idler photons. Since a photon's energy is inversely
    proportional to its wavelength (E = hc/λ), this relationship can be
    expressed in terms of wavelengths:

    1/λ_pump = 1/λ_signal + 1/λ_idler

    This script calculates the idler wavelength given the pump and signal
    wavelengths.

    Note: The question of whether a material *can* exhibit SPDC depends on its
    crystal symmetry. Materials with inversion symmetry (centrosymmetric),
    like the common phases of free-standing boron nanosheets, have a vanishing
    second-order nonlinearity and are not expected to exhibit this effect.
    This script simply demonstrates the math if the effect were to occur.
    """

    # --- User-defined values ---
    # Example values are in nanometers (nm)
    pump_wavelength = 405  # A common laser wavelength (blue/violet)
    signal_wavelength = 810 # A common case is degenerate SPDC where signal/idler are double the pump wavelength

    # --- Calculation ---
    # The equation is: 1/λ_idler = 1/λ_pump - 1/λ_signal
    # First, check if the signal wavelength is physically possible (must be longer than the pump)
    if signal_wavelength <= pump_wavelength:
        print(f"Error: Signal wavelength ({signal_wavelength} nm) must be greater than pump wavelength ({pump_wavelength} nm).")
        return

    # Calculate the idler wavelength
    idler_wavelength = 1 / (1 / pump_wavelength - 1 / signal_wavelength)

    # --- Output ---
    print("Illustrating Energy Conservation in Spontaneous Parametric Down-Conversion (SPDC):")
    print("-" * 75)
    print(f"Given a pump photon wavelength of {pump_wavelength} nm and a signal photon wavelength of {signal_wavelength} nm...")
    print("\nWe use the conservation of energy formula: ")
    print("1/λ_pump = 1/λ_signal + 1/λ_idler")
    
    print("\nPlugging in the numbers to solve for the idler wavelength:")
    print(f"1/{pump_wavelength} = 1/{signal_wavelength} + 1/λ_idler")
    
    # Rearranging the formula to solve for the idler
    print(f"1/λ_idler = 1/{pump_wavelength} - 1/{signal_wavelength}")
    
    print("\nResult:")
    # Using round to present a cleaner number
    print(f"The calculated idler photon wavelength is: {round(idler_wavelength, 2)} nm")
    
    print("\nFinal Verifying Equation:")
    # Showing the full equation with the calculated idler value
    print(f"1/{pump_wavelength} = 1/{signal_wavelength} + 1/{round(idler_wavelength, 2)}")


# Execute the function
calculate_spdc_wavelengths()