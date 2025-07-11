def check_visibility():
    """
    This function checks if the light from Mercury's sodium tail would be visible
    through a specified bandpass filter.
    """
    # Define the wavelengths in nanometers (nm)
    sodium_emission_wavelength = 589.3  # Average wavelength of the sodium D-lines, which gives the tail its color
    filter_center = 486               # The center wavelength of the filter used
    filter_bandwidth = 1              # The width of the filter's passband

    # Calculate the filter's effective range
    filter_min = filter_center - (filter_bandwidth / 2)
    filter_max = filter_center + (filter_bandwidth / 2)

    print("Step 1: Define emission and filter properties.")
    print(f"Emission wavelength of Sodium tail: {sodium_emission_wavelength} nm")
    print(f"Filter is centered at: {filter_center} nm")
    print(f"Filter bandwidth is: {filter_bandwidth} nm")
    
    print("\nStep 2: Calculate the filter's passband range.")
    print(f"The filter allows light from {filter_min} nm to {filter_max} nm to pass through.")

    print("\nStep 3: Check for visibility.")
    # The 'equation' is the boolean check below
    is_visible = filter_min <= sodium_emission_wavelength <= filter_max
    
    print(f"Is {filter_min} <= {sodium_emission_wavelength} <= {filter_max}? The answer is {is_visible}.")

    if is_visible:
        print("\nConclusion: The sodium tail is VISIBLE through the filter.")
    else:
        print(f"\nConclusion: The light from the sodium tail ({sodium_emission_wavelength} nm) is outside the filter's passband ({filter_min}-{filter_max} nm).")
        print("Therefore, no tail would be visible.")

# Run the function
check_visibility()