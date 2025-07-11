def analyze_observation():
    """
    Analyzes the visibility of Mercury's sodium tail based on the provided instrument details.
    """
    # Define the key wavelengths in nanometers (nm)
    sodium_emission_line = 589.0  # The primary emission wavelength of the sodium tail (yellow)
    filter_center_wavelength = 486  # The center of the instrument's filter (blue)
    filter_bandwidth = 1.0         # The width of the filter

    # Calculate the filter's operational range
    filter_min_pass = filter_center_wavelength - (filter_bandwidth / 2)
    filter_max_pass = filter_center_wavelength + (filter_bandwidth / 2)

    # Print the parameters of the problem for clarity
    print(f"Step 1: Define the light source.")
    print(f"The sodium tail's characteristic emission is at ~{sodium_emission_line} nm (yellow).")
    print("-" * 30)

    print(f"Step 2: Define the instrument's filter.")
    print(f"The filter is centered at {filter_center_wavelength} nm (blue).")
    print(f"The filter allows light from {filter_min_pass} nm to {filter_max_pass} nm to pass through.")
    print("-" * 30)

    # Step 3: Check if the sodium light can pass through the filter
    print("Step 3: Compare the light source to the filter.")
    can_be_seen = filter_min_pass <= sodium_emission_line <= filter_max_pass

    if can_be_seen:
        print(f"Result: The sodium emission at {sodium_emission_line} nm falls within the filter's range.")
        print("Conclusion: The tail would be visible.")
    else:
        print(f"Result: The sodium emission at {sodium_emission_line} nm is outside the filter's range of {filter_min_pass}-{filter_max_pass} nm.")
        print("\nConclusion: The filter blocks the light from the sodium tail. Therefore, no tail would be visible, regardless of Mercury's position (aphelion or perihelion).")

analyze_observation()