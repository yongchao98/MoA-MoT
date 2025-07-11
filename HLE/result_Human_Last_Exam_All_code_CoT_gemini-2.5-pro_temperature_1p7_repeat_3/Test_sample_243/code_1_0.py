def check_visibility():
    """
    This function analyzes the visibility of Mercury's sodium tail based on observation parameters.
    """
    # Wavelength of the light emitted by Mercury's sodium tail in nanometers
    sodium_emission_wavelength = 589

    # Wavelength the observation filter is centered on in nanometers
    filter_wavelength = 486

    print(f"The filter is looking for light at a wavelength of {filter_wavelength} nm.")
    print(f"Mercury's sodium tail emits light at a wavelength of approximately {sodium_emission_wavelength} nm.")

    # Check if the emission is within the filter's passband.
    # For this problem, we assume the bandpass is very narrow and doesn't include the emission.
    if filter_wavelength != sodium_emission_wavelength:
        print("\nConclusion: The filter blocks the light from the sodium tail.")
        print("Therefore, no tail would be visible.")
    else:
        # This case is not met in the problem.
        print("The tail would be visible.")

check_visibility()