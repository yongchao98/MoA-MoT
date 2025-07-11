def check_visibility():
    """
    Checks if Mercury's sodium tail is visible through a specific filter.
    """
    # Wavelength of the prominent sodium D-line emission from Mercury's tail (in nm)
    sodium_emission_wavelength = 589.0

    # Properties of the bandpass filter used for observation (in nm)
    filter_center_wavelength = 486.0
    filter_width = 1.0

    # Calculate the range of wavelengths the filter allows to pass
    filter_min = filter_center_wavelength - filter_width / 2.0
    filter_max = filter_center_wavelength + filter_width / 2.0

    print(f"The sodium tail emits light at approximately {sodium_emission_wavelength} nm (yellow).")
    print(f"The telescope's filter is centered at {filter_center_wavelength} nm with a {filter_width} nm width.")
    print(f"This filter only allows light between {filter_min} nm and {filter_max} nm to be seen.")

    # Check if the sodium emission is within the filter's range
    if filter_min <= sodium_emission_wavelength <= filter_max:
        print("\nConclusion: The sodium emission is within the filter's range. The tail would be visible.")
    else:
        print(f"\nConclusion: The sodium emission at {sodium_emission_wavelength} nm is outside the filter's range.")
        print("Therefore, the filter blocks the light from the sodium tail, and no tail would be visible.")

check_visibility()