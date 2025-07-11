def check_visibility():
    """
    Determines if Mercury's sodium tail is visible through a specific filter.
    """
    # Wavelength of the prominent sodium D-line emission from the tail
    sodium_emission_wavelength = 589.0  # in nanometers (nm), which is yellow light

    # Properties of the telescope's filter
    filter_center_wavelength = 486.0  # in nm, which is blue-green light
    filter_bandwidth = 1.0            # in nm

    # Calculate the wavelength range the filter allows to pass
    filter_min_wavelength = filter_center_wavelength - (filter_bandwidth / 2)
    filter_max_wavelength = filter_center_wavelength + (filter_bandwidth / 2)

    print(f"The sodium tail emits light primarily at a wavelength of {sodium_emission_wavelength} nm.")
    print(f"The telescope's filter is centered at {filter_center_wavelength} nm and has a bandwidth of {filter_bandwidth} nm.")
    print(f"This means the filter only allows light from {filter_min_wavelength} nm to {filter_max_wavelength} nm to pass through.")

    # Check if the sodium emission falls within the filter's range
    if sodium_emission_wavelength >= filter_min_wavelength and sodium_emission_wavelength <= filter_max_wavelength:
        print("\nResult: The tail would be visible.")
    else:
        print(f"\nResult: The filter blocks the {sodium_emission_wavelength} nm light from the sodium tail. Therefore, no tail would be visible.")

check_visibility()