import sys

def analyze_visibility():
    """
    Analyzes the visibility of Mercury's sodium tail given specific
    observational parameters.
    """
    # Define the physical and observational parameters
    sodium_emission_wavelength = 589.0  # Wavelength of sodium D-lines in nm (yellow)
    filter_center = 486.0               # Center of the bandpass filter in nm (blue-green)
    filter_width = 1.0                  # Width of the bandpass filter in nm

    # Calculate the filter's passband (the range of wavelengths it allows through)
    filter_lower_bound = filter_center - (filter_width / 2.0)
    filter_upper_bound = filter_center + (filter_width / 2.0)

    # Print the analysis steps
    print(f"Step 1: The sodium tail emits light primarily at {sodium_emission_wavelength} nm.")
    print(f"Step 2: The telescope uses a filter that only allows light between {filter_lower_bound} nm and {filter_upper_bound} nm to pass through.")

    # Check if the sodium emission falls within the filter's passband
    is_visible = filter_lower_bound <= sodium_emission_wavelength <= filter_upper_bound
    
    # Use an "equation" format for the comparison check as requested
    print(f"\nStep 3: Check if the emission is in the filter's range.")
    print(f"Is {filter_lower_bound} <= {sodium_emission_wavelength} <= {filter_upper_bound}?")
    print(f"Result of check: {is_visible}")


    # Final conclusion based on the check
    if not is_visible:
        print("\nConclusion: The light from the sodium tail ({sodium_emission_wavelength} nm) is blocked by the filter ({filter_center} nm). Therefore, no tail would be visible.")
    else:
        # This case won't be reached with the given numbers
        print("\nConclusion: The light from the sodium tail would be visible through the filter.")

if __name__ == '__main__':
    analyze_visibility()