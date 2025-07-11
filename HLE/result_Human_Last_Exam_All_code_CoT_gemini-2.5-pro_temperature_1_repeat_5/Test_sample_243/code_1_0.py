def check_visibility():
    """
    This function determines if Mercury's sodium tail would be visible
    through a specific bandpass filter.
    """
    # Wavelength of the primary sodium D-line emission from the tail (in nanometers)
    sodium_emission_wavelength = 589.0

    # Properties of the bandpass filter (in nanometers)
    filter_center_wavelength = 486.0
    filter_bandwidth = 1.0

    # Calculate the filter's transmission range
    filter_min_wavelength = filter_center_wavelength - (filter_bandwidth / 2)
    filter_max_wavelength = filter_center_wavelength + (filter_bandwidth / 2)

    print(f"Step 1: The tail is composed of sodium, which emits light primarily at {sodium_emission_wavelength} nm (yellow).")
    print(f"Step 2: The telescope uses a filter centered at {filter_center_wavelength} nm with a {filter_bandwidth} nm bandwidth.")
    print(f"This means the filter only passes light between {filter_min_wavelength} nm and {filter_max_wavelength} nm.")

    # Check if the sodium emission is within the filter's range
    if filter_min_wavelength <= sodium_emission_wavelength <= filter_max_wavelength:
        is_visible = True
    else:
        is_visible = False

    print("\n--- Calculation ---")
    print(f"Is {sodium_emission_wavelength} nm between {filter_min_wavelength} nm and {filter_max_wavelength} nm? {is_visible}")

    print("\n--- Conclusion ---")
    if not is_visible:
        print("The yellow light from the sodium tail is outside the filter's range and will be blocked.")
        print("Therefore, no tail would be visible, regardless of its length.")
    else:
        # This case won't be reached with the given numbers
        print("The light from the sodium tail would be visible.")

# Run the analysis
check_visibility()
<<<B>>>