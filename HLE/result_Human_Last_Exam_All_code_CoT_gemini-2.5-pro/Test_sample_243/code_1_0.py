def analyze_observation():
    """
    Analyzes the visibility of Mercury's sodium tail with a specific filter.
    """
    # Wavelength of the filter in nanometers (nm)
    filter_center_wavelength = 486
    # Width of the filter's bandpass in nm
    filter_width = 1

    # Wavelength of the primary sodium emission lines (D-lines) in nm
    sodium_emission_wavelength = 589

    # Calculate the filter's passband range
    filter_passband_min = filter_center_wavelength - (filter_width / 2.0)
    filter_passband_max = filter_center_wavelength + (filter_width / 2.0)

    # Print the parameters of the observation
    print(f"Observation Parameters:")
    print(f"  - The filter is centered at {filter_center_wavelength} nm (blue-cyan light).")
    print(f"  - The filter's passband is from {filter_passband_min} nm to {filter_passband_max} nm.")
    print(f"  - Mercury's tail is composed of sodium, which emits strongly at ~{sodium_emission_wavelength} nm (yellow light).")
    print("-" * 30)

    # Check if the sodium emission is visible through the filter
    if sodium_emission_wavelength >= filter_passband_min and sodium_emission_wavelength <= filter_passband_max:
        print("Conclusion: The sodium tail would be visible.")
    else:
        print("Conclusion: The sodium emission wavelength is outside the filter's passband.")
        print("Therefore, no tail would be visible.")

# Run the analysis
analyze_observation()
