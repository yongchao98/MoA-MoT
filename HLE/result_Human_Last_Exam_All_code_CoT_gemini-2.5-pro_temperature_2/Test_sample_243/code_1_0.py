import sys

def check_visibility():
    """
    This function checks if the light from Mercury's sodium tail is visible
    through the specified optical filter.
    """
    # Wavelength of the light emitted by the sodium tail (Sodium D-lines)
    sodium_wavelength_nm = 589.0

    # Properties of the bandpass filter
    filter_center_nm = 486.0
    filter_width_nm = 1.0

    # Calculate the range of the filter
    filter_min_nm = filter_center_nm - (filter_width_nm / 2.0)
    filter_max_nm = filter_center_nm + (filter_width_nm / 2.0)

    # Print the analysis steps with the numbers involved
    print("Analysis of Observation:")
    print("-" * 25)
    print(f"1. The primary emission wavelength from the sodium tail is: {sodium_wavelength_nm} nm")
    print(f"2. The observer is using a filter centered at {filter_center_nm} nm with a width of {filter_width_nm} nm.")
    print(f"3. This filter's equation for its passband is: {filter_min_nm} nm <= Wavelength <= {filter_max_nm} nm")

    # Final check: Does the sodium wavelength fall within the filter's range?
    if filter_min_nm <= sodium_wavelength_nm <= filter_max_nm:
        print("\nConclusion: The sodium light falls within the filter's range and would be visible.")
    else:
        print(f"\nConclusion: The sodium light at {sodium_wavelength_nm} nm is outside the filter's range.")
        print("Therefore, no tail would be visible.")


if __name__ == '__main__':
    # Some platforms might have issues with unicode characters.
    # We will try to set the encoding to UTF-8 to prevent errors.
    if sys.stdout.encoding != 'utf-8':
        try:
            sys.stdout.reconfigure(encoding='utf-8')
        except TypeError:
            # In some environments, reconfigure might not be available.
            # We will ignore the error and proceed.
            pass
    check_visibility()