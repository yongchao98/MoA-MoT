def check_visibility():
    """
    Checks if Mercury's sodium tail is visible through a specific filter.
    """
    # Wavelength of the prominent sodium D-lines (in nanometers)
    sodium_d_line_wavelength = 589.0

    # Filter specifications (in nanometers)
    filter_center_wavelength = 486.0
    filter_bandwidth = 1.0

    # Calculate the filter's passband
    min_wavelength = filter_center_wavelength - (filter_bandwidth / 2.0)
    max_wavelength = filter_center_wavelength + (filter_bandwidth / 2.0)

    print(f"The primary emission wavelength of Mercury's sodium tail is {sodium_d_line_wavelength} nm (yellow light).")
    print(f"The telescope's filter allows light between {min_wavelength} nm and {max_wavelength} nm to pass through.")

    # Check if the sodium emission is within the filter's passband
    if min_wavelength <= sodium_d_line_wavelength <= max_wavelength:
        print("\nConclusion: The sodium emission line is within the filter's passband.")
        print("The sodium tail would be visible.")
    else:
        print(f"\nConclusion: The sodium emission line at {sodium_d_line_wavelength} nm is NOT within the filter's passband of {min_wavelength}-{max_wavelength} nm.")
        print("Therefore, no significant light from the sodium tail can reach the detector.")
        print("The sodium tail would not be visible with this instrumental setup.")
        print("\nThis corresponds to answer choice B.")

check_visibility()