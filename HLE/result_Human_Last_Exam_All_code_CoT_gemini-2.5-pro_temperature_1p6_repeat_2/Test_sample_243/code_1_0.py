import sys

def analyze_mercury_observation():
    """
    Analyzes the visibility of Mercury's sodium tail based on
    the emission wavelength and the filter used for observation.
    """
    # Define the physical constants for the problem
    # Wavelength of Mercury's sodium tail emission (D-lines) in nanometers (nm)
    sodium_emission_wavelength_nm = 589.0
    
    # Wavelength of the center of the bandpass filter in nanometers (nm)
    filter_center_wavelength_nm = 486.0
    
    # Width of the bandpass filter in nanometers (nm)
    filter_bandwidth_nm = 1.0

    # Calculate the filter's effective range
    filter_min_pass_nm = filter_center_wavelength_nm - (filter_bandwidth_nm / 2.0)
    filter_max_pass_nm = filter_center_wavelength_nm + (filter_bandwidth_nm / 2.0)

    # Print the analysis
    print("--- Analysis of Mercury's Tail Observation ---")
    print(f"1. The primary light from Mercury's tail is sodium emission, which occurs at a wavelength of {sodium_emission_wavelength_nm} nm (yellow light).")
    print(f"2. The observation is made with a filter centered at {filter_center_wavelength_nm} nm (blue-green light) with a bandwidth of {filter_bandwidth_nm} nm.")
    
    # We construct the equation to check if the emission is in the filter's range.
    # The equation is: filter_min <= emission <= filter_max
    # Substituting the numbers: filter_min_pass_nm <= sodium_emission_wavelength_nm <= filter_max_pass_nm
    print("\n3. To be visible, the sodium emission must fall within the filter's passband.")
    print(f"   Checking if {sodium_emission_wavelength_nm} nm is between {filter_min_pass_nm} nm and {filter_max_pass_nm} nm.")
    
    # Perform the logical check
    if sodium_emission_wavelength_nm >= filter_min_pass_nm and sodium_emission_wavelength_nm <= filter_max_pass_nm:
        # This case is physically incorrect for this problem but included for completeness
        print("\nConclusion: The sodium emission is within the filter's range. The tail would be visible.")
    else:
        # The correct conclusion based on the numbers
        print("\nConclusion: The sodium emission is NOT within the filter's range.")
        print("The filter, which passes blue-green light, completely blocks the yellow light from the sodium tail.")
        print("Therefore, no tail would be visible in the image.")

# Run the analysis
analyze_mercury_observation()