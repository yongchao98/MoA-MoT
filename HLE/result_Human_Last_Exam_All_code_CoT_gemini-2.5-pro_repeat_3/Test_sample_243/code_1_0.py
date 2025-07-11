def analyze_observation():
    """
    Analyzes the visibility of Mercury's sodium tail based on the filter used.
    """
    # Define the primary emission wavelength of the sodium tail in nanometers (nm)
    sodium_emission_wavelength = 589.0

    # Define the properties of the optical filter being used in nm
    filter_center_wavelength = 486.0
    filter_bandpass_width = 1.0

    # Calculate the range of wavelengths the filter allows to pass
    filter_min_wavelength = filter_center_wavelength - (filter_bandpass_width / 2.0)
    filter_max_wavelength = filter_center_wavelength + (filter_bandpass_width / 2.0)

    print("Step 1: Identify the light source.")
    print(f"Mercury's tail is made of sodium, which emits light primarily at {sodium_emission_wavelength} nm (yellow).")
    print("-" * 50)
    
    print("Step 2: Identify the observation filter's properties.")
    print(f"The telescope is using a filter centered at {filter_center_wavelength} nm (blue-green).")
    print(f"The filter has a bandpass width of {filter_bandpass_width} nm.")
    print(f"This means the filter only allows light from {filter_min_wavelength} nm to {filter_max_wavelength} nm to pass through.")
    print("-" * 50)

    print("Step 3: Compare the light source to the filter.")
    # Check if the sodium emission falls within the filter's range
    if sodium_emission_wavelength >= filter_min_wavelength and sodium_emission_wavelength <= filter_max_wavelength:
        is_visible = True
    else:
        is_visible = False
    
    print(f"Is the sodium emission wavelength ({sodium_emission_wavelength} nm) within the filter's range ({filter_min_wavelength}-{filter_max_wavelength} nm)? {is_visible}")
    print("-" * 50)

    print("Final Conclusion:")
    if not is_visible:
        print("The filter blocks the light from the sodium tail.")
        print("Therefore, no tail would be visible, regardless of its actual length.")
    else:
        print("The tail would be visible through the filter.")

analyze_observation()