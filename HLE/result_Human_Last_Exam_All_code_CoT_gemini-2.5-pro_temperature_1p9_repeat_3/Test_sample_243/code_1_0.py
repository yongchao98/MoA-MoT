def mercury_tail_observation():
    """
    Analyzes the visibility of Mercury's sodium tail with a specific filter.
    """
    
    # Wavelength of the prominent sodium D-lines emission from Mercury's tail
    sodium_emission_wavelength_nm = 589.0 
    
    # Color corresponding to the sodium emission wavelength
    sodium_emission_color = "yellow"
    
    # The central wavelength of the bandpass filter used for observation
    filter_wavelength_nm = 486.0
    
    # The color corresponding to the filter's wavelength (Hydrogen-beta line)
    filter_color = "blue-green"
    
    # The width of the filter's bandpass
    filter_width_nm = 1.0
    
    # Calculate the filter's effective range
    min_pass_wavelength = filter_wavelength_nm - filter_width_nm / 2
    max_pass_wavelength = filter_wavelength_nm + filter_width_nm / 2

    print("Step 1: The tail is made of sodium, which glows with a characteristic color.")
    print(f"         - Sodium tail emission wavelength: {sodium_emission_wavelength_nm} nm ({sodium_emission_color}).")
    
    print("\nStep 2: The observation is made using a specific filter.")
    print(f"         - The filter is centered at: {filter_wavelength_nm} nm ({filter_color}).")
    print(f"         - The filter allows light between {min_pass_wavelength} nm and {max_pass_wavelength} nm to pass through.")

    print("\nStep 3: Compare the tail's emission to the filter's passband.")
    
    # Check if the sodium emission falls within the filter's range
    if min_pass_wavelength <= sodium_emission_wavelength_nm <= max_pass_wavelength:
        print("         - The filter allows the sodium light to pass. The tail would be visible.")
    else:
        print("         - The filter blocks the sodium light.")
        print("         - The emission wavelength (589.0 nm) is outside the filter's range (485.5-486.5 nm).")

    print("\nConclusion: Because the filter blocks the light emitted by the sodium tail, the tail would not be visible, regardless of Mercury's position in its orbit.")

mercury_tail_observation()