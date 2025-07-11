def explain_mercury_observation():
    """
    Explains the visibility of Mercury's sodium tail with a 486nm filter.
    """
    # Wavelength of the prominent Sodium D-line emissions from Mercury's tail
    sodium_emission_wavelength_nm = 589.0

    # The central wavelength of the filter being used
    filter_center_wavelength_nm = 486.0

    # The width of the filter's passband
    filter_width_nm = 1.0

    # Calculate the filter's range
    filter_min = filter_center_wavelength_nm - filter_width_nm / 2
    filter_max = filter_center_wavelength_nm + filter_width_nm / 2

    print("Step 1: Understand the source of light from Mercury's sodium tail.")
    print(f"The sodium tail shines primarily because sodium atoms emit light at specific wavelengths.")
    print(f"The most significant of these wavelengths are the Sodium D-lines at ~{sodium_emission_wavelength_nm} nm.")
    print("-" * 30)

    print("Step 2: Understand the instrument being used.")
    print(f"The observation uses a narrow bandpass filter centered at {filter_center_wavelength_nm} nm.")
    print(f"This filter only allows light between {filter_min} nm and {filter_max} nm to pass through.")
    print("-" * 30)

    print("Step 3: Compare the light source and the filter.")
    print(f"The sodium tail's light is at {sodium_emission_wavelength_nm} nm.")
    print(f"The filter's range is {filter_min}-{filter_max} nm.")

    # Check if the sodium emission is within the filter's range
    if filter_min <= sodium_emission_wavelength_nm <= filter_max:
        print("Conclusion: The filter would pass the light from the sodium tail.")
    else:
        print("Conclusion: The filter will block the light from the sodium tail because its wavelength is outside the filter's range.")
        print("Therefore, the sodium tail would not be visible at all, regardless of Mercury's position (aphelion or perihelion).")

# Execute the explanation
explain_mercury_observation()