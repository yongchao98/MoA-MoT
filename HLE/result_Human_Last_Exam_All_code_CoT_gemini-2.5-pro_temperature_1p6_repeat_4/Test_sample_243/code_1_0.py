def check_visibility():
    """
    Analyzes whether Mercury's sodium tail would be visible given the
    instrument parameters.
    """
    # Parameters from the problem description
    filter_center_wl = 486  # in nanometers
    filter_width = 1        # in nanometers
    sodium_emission_wl = 589.0  # Approximate wavelength of sodium D-lines in nm

    # Calculate the wavelength range the filter allows through
    min_pass_wl = filter_center_wl - (filter_width / 2)
    max_pass_wl = filter_center_wl + (filter_width / 2)

    # Output the parameters for clarity
    print(f"Observation Filter Center: {filter_center_wl} nm")
    print(f"Filter Passband: {min_pass_wl} nm to {max_pass_wl} nm")
    print(f"Sodium Tail Emission Wavelength: ~{sodium_emission_wl} nm")
    
    # Check if the sodium emission is within the filter's range
    if sodium_emission_wl >= min_pass_wl and sodium_emission_wl <= max_pass_wl:
        print("\nResult: The sodium emission is within the filter's range.")
    else:
        print("\nResult: The filter blocks the light from the sodium tail.")
        print("Therefore, no tail would be visible.")

check_visibility()