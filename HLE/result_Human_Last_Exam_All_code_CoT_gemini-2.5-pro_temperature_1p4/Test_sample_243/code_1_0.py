import sys

def analyze_observation():
    """
    Analyzes the visibility of Mercury's sodium tail based on emission and filter wavelengths.
    """
    # Key physical parameters in nanometers (nm)
    sodium_emission_wavelength = 589
    filter_center_wavelength = 486
    filter_bandwidth = 1

    # The problem describes a comparison between the light emitted by the tail
    # and the light allowed through the filter.
    # The tail's visibility depends on whether its emission wavelength falls
    # within the filter's passband.
    
    print("--- Analysis of Mercury's Tail Observation ---")
    print(f"1. The sodium tail emits light primarily at the sodium D-line wavelength.")
    print(f"   Sodium Emission Wavelength = {sodium_emission_wavelength} nm (Yellow)")
    
    print("\n2. The telescope uses a specific bandpass filter.")
    print(f"   Filter Center Wavelength = {filter_center_wavelength} nm (Blue-Green)")
    print(f"   Filter Bandwidth = {filter_bandwidth} nm")

    # Calculate the filter's operational range
    filter_min = filter_center_wavelength - (filter_bandwidth / 2)
    filter_max = filter_center_wavelength + (filter_bandwidth / 2)

    print(f"\n3. This means the filter only allows light between {filter_min} nm and {filter_max} nm to pass through.")

    # Check for overlap
    print("\n--- Conclusion ---")
    if sodium_emission_wavelength >= filter_min and sodium_emission_wavelength <= filter_max:
        print(f"The sodium emission at {sodium_emission_wavelength} nm IS within the filter's range.")
        print("Therefore, the tail would be visible.")
    else:
        print(f"The sodium emission at {sodium_emission_wavelength} nm IS NOT within the filter's range of {filter_min}-{filter_max} nm.")
        print("Therefore, the light from the sodium tail will be blocked and no tail would be visible.")

if __name__ == '__main__':
    analyze_observation()
