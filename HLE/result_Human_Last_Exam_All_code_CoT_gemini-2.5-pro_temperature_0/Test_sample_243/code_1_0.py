def analyze_mercury_observation():
    """
    Analyzes the visibility of Mercury's sodium tail based on observation parameters.
    """
    # Wavelength of the prominent sodium D-lines emission from Mercury's tail
    sodium_d1_wavelength_nm = 589.0
    sodium_d2_wavelength_nm = 589.6

    # Parameters of the bandpass filter used for observation
    filter_center_nm = 486
    filter_width_nm = 1

    # Calculate the range of the filter
    filter_min_nm = filter_center_nm - (filter_width_nm / 2)
    filter_max_nm = filter_center_nm + (filter_width_nm / 2)

    print("Analysis of Mercury's Sodium Tail Observation:")
    print("-" * 50)
    print(f"1. Mercury's sodium tail primarily emits yellow light at the sodium D-lines.")
    print(f"   - Sodium D1 line: {sodium_d1_wavelength_nm} nm")
    print(f"   - Sodium D2 line: {sodium_d2_wavelength_nm} nm")
    print("\n2. The observation is conducted using a specific bandpass filter.")
    print(f"   - Filter Center: {filter_center_nm} nm")
    print(f"   - Filter Width: {filter_width_nm} nm")
    print(f"   - This filter only allows light between {filter_min_nm} nm and {filter_max_nm} nm to pass through.")
    print("\n3. Comparing the emission and the filter:")

    # Check if the sodium light can pass through the filter
    is_visible = (filter_min_nm <= sodium_d1_wavelength_nm <= filter_max_nm) or \
                 (filter_min_nm <= sodium_d2_wavelength_nm <= filter_max_nm)

    if not is_visible:
        print(f"   - The sodium emission at ~{int((sodium_d1_wavelength_nm + sodium_d2_wavelength_nm)/2)} nm is completely outside the filter's range of {filter_min_nm}-{filter_max_nm} nm.")
        print("\nConclusion:")
        print("The light from the sodium tail will be blocked by the filter.")
        print("Therefore, regardless of Mercury's position (aphelion or perihelion), no tail would be visible with this instrument setup.")
    else:
        # This case will not be reached with the given numbers
        print("The sodium tail would be visible.")

analyze_mercury_observation()