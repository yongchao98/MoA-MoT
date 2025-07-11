import math

def check_visibility():
    """
    Analyzes the visibility of Mercury's sodium tail based on the provided filter.
    """
    # Wavelength of the prominent sodium D-line emissions (in nanometers)
    sodium_emission_wavelength = 589.0
    sodium_emission_color = "yellow"

    # Properties of the bandpass filter (in nanometers)
    filter_center_wavelength = 486.0
    filter_width = 1.0
    filter_color = "blue-green (H-beta)"

    # Calculate the range of the filter
    filter_min = filter_center_wavelength - (filter_width / 2.0)
    filter_max = filter_center_wavelength + (filter_width / 2.0)

    print("--- Analysis of Observation ---")
    print(f"Mercury's sodium tail emits light primarily at ~{sodium_emission_wavelength}nm, which appears {sodium_emission_color}.")
    print(f"The telescope filter is centered at {filter_center_wavelength}nm and has a passband from {filter_min}nm to {filter_max}nm.")

    # Check if the sodium emission is within the filter's range
    if sodium_emission_wavelength >= filter_min and sodium_emission_wavelength <= filter_max:
        print("\nConclusion: The sodium emission IS WITHIN the filter's passband. The tail would be visible.")
    else:
        print(f"\nConclusion: The sodium emission at {sodium_emission_wavelength}nm IS NOT WITHIN the filter's passband of {filter_min}-{filter_max}nm.")
        print("The filter will block the light from the sodium tail.")
        print("Therefore, no tail would be visible.")
        print("\nCorrect Answer Choice: B")

if __name__ == '__main__':
    check_visibility()
