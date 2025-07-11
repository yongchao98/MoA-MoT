def check_visibility():
    """
    Analyzes the visibility of Mercury's sodium tail based on filter specifications.
    """
    # Wavelength of the prominent sodium D-lines (in nanometers)
    sodium_emission_wavelength_1 = 589.0
    sodium_emission_wavelength_2 = 589.6
    
    # Filter specifications (in nanometers)
    filter_center_wavelength = 486.0
    filter_bandpass_width = 1.0
    
    # Calculate the range of wavelengths the filter allows to pass
    filter_min_wavelength = filter_center_wavelength - (filter_bandpass_width / 2)
    filter_max_wavelength = filter_center_wavelength + (filter_bandpass_width / 2)
    
    print("Analysis of Mercury's Sodium Tail Visibility:")
    print("-" * 45)
    print(f"The sodium tail primarily emits light at {sodium_emission_wavelength_1} nm and {sodium_emission_wavelength_2} nm (yellow).")
    print(f"The telescope filter is centered at {filter_center_wavelength} nm (blue-green).")
    print(f"The filter allows light in the range {filter_min_wavelength:.1f} nm to {filter_max_wavelength:.1f} nm to pass through.")
    
    # Check if the sodium emission is within the filter's range
    is_visible_1 = filter_min_wavelength <= sodium_emission_wavelength_1 <= filter_max_wavelength
    is_visible_2 = filter_min_wavelength <= sodium_emission_wavelength_2 <= filter_max_wavelength
    
    print("\nConclusion:")
    if not is_visible_1 and not is_visible_2:
        print("The sodium emission wavelengths are outside the filter's passband.")
        print("Therefore, the light from the sodium tail will be blocked by the filter.")
        print("Result: No tail would be visible.")
    else:
        print("The sodium emission is within the filter's passband.")
        print("Result: The tail would be visible.")

check_visibility()