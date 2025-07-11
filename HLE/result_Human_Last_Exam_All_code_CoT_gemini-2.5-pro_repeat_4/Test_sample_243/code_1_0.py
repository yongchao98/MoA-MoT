import sys

def analyze_mercury_observation():
    """
    Analyzes the visibility of Mercury's sodium tail based on the provided
    observational filter parameters.
    """
    # Key physical and observational parameters
    sodium_emission_wavelength_nm = 589.0  # Main emission line for sodium (in nanometers), which is yellow.
    filter_center_wavelength_nm = 486.0    # Filter centered at H-beta line (in nanometers), which is blue.
    filter_bandwidth_nm = 1.0              # The width of the filter passband.

    # Calculate the range of wavelengths the filter allows to pass
    filter_min_nm = filter_center_wavelength_nm - (filter_bandwidth_nm / 2)
    filter_max_nm = filter_center_wavelength_nm + (filter_bandwidth_nm / 2)

    print("--- Analysis of Mercury's Sodium Tail Observation ---")
    print(f"1. Mercury's sodium tail primarily emits light at the sodium D-line, which is around {sodium_emission_wavelength_nm} nm (yellow).")
    print(f"2. The observation is being made with a filter centered at {filter_center_wavelength_nm} nm.")
    print(f"3. With a bandwidth of {filter_bandwidth_nm} nm, this filter only allows light between {filter_min_nm} nm and {filter_max_nm} nm to pass through.")

    # Check if the sodium emission is within the filter's passband
    if filter_min_nm <= sodium_emission_wavelength_nm <= filter_max_nm:
        # This case will not be reached with the given numbers
        print("\nResult: The sodium emission falls within the filter's range. The tail would be visible.")
    else:
        print(f"\nResult: The light from the sodium tail ({sodium_emission_wavelength_nm} nm) is completely outside the filter's passband ({filter_min_nm}-{filter_max_nm} nm).")
        print("The filter will block the light from the sodium tail.")

    print("\n--- Conclusion ---")
    print("Regardless of whether Mercury is at aphelion or perihelion, the filter being used is incorrect for observing the sodium tail.")
    print("Therefore, no tail would be visible.")

# Execute the analysis
analyze_mercury_observation()
