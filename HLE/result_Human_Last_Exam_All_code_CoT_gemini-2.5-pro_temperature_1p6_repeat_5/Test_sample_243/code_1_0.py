import sys

def check_visibility():
    """
    Analyzes the visibility of Mercury's sodium tail through a specific filter.
    """
    # Wavelength of the prominent sodium D-lines (in nanometers)
    sodium_emission_wavelength_nm = 589.0

    # Properties of the bandpass filter
    filter_center_nm = 486
    filter_width_nm = 1

    # Calculate the filter's passband range
    filter_lower_bound = filter_center_nm - (filter_width_nm / 2)
    filter_upper_bound = filter_center_nm + (filter_width_nm / 2)

    print("Analyzing the visibility of Mercury's sodium tail:")
    print(f"1. The sodium tail primarily emits yellow light at ~{sodium_emission_wavelength_nm:.1f} nm.")
    print(f"2. The telescope filter only allows light from {filter_lower_bound:.1f} nm to {filter_upper_bound:.1f} nm to pass through (centered at {filter_center_nm} nm).")

    # Check if the sodium emission is within the filter's range
    if filter_lower_bound <= sodium_emission_wavelength_nm <= filter_upper_bound:
        print("\nResult: The sodium emission wavelength IS WITHIN the filter's range. The tail would be visible.")
    else:
        print(f"\nResult: The sodium emission at {sodium_emission_wavelength_nm:.1f} nm is NOT WITHIN the filter's range of {filter_lower_bound:.1f}-{filter_upper_bound:.1f} nm.")
        print("Therefore, the filter will block the light from the sodium tail, and no tail would be visible.")
    
    print("\nAdditionally, while the tail's length varies (it is shorter at aphelion and longer at perihelion), this is irrelevant because it cannot be seen with this filter.")

# The problem is conceptual, but the code below can run to show the logic.
# To prevent execution in certain environments, it's wrapped in a main guard.
if __name__ == '__main__':
    # Running this script prints the reasoning. The final choice is based on this logic.
    if len(sys.argv) > 1 and sys.argv[1] == '--execute':
        check_visibility()

# Based on the analysis that the filter blocks the light, the correct answer is that no tail would be visible.
