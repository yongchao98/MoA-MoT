import sys

def analyze_mercury_observation():
    """
    Analyzes the visibility of Mercury's sodium tail based on observation parameters.
    """
    # Step 1: Define the observational parameters from the problem
    filter_center_wl_nm = 486
    filter_bandpass_nm = 1
    
    # Step 2: Define the physical properties of Mercury's sodium tail
    # The tail is primarily made of sodium, which has strong emission lines (D-lines)
    sodium_emission_wl_nm = 589.0 # We'll use the stronger D2 line for this check
    
    # Step 3: Calculate the wavelength range the filter allows to pass
    filter_min_wl = filter_center_wl_nm - (filter_bandpass_nm / 2)
    filter_max_wl = filter_center_wl_nm + (filter_bandpass_nm / 2)

    # Step 4: Systematically print the reasoning
    print("--- Analysis of the Observation ---")
    print(f"1. The observation uses a filter centered at {filter_center_wl_nm} nm with a {filter_bandpass_nm} nm bandpass.")
    # The final print needs each number in the final equation. So we are showing how the range is calculated.
    print(f"   Calculation for the filter's pass-through range: {filter_center_wl_nm} nm +/- ({filter_bandpass_nm}/2) nm")
    print(f"   This means the filter only lets light between {filter_min_wl} nm and {filter_max_wl} nm pass through.")
    print("\n2. Mercury's tail is composed of sodium atoms, which strongly emit yellow light.")
    print(f"   The primary sodium emission wavelength is ~{sodium_emission_wl_nm} nm.")
    
    # Step 5: Check if the sodium emission is within the filter's range
    is_visible = filter_min_wl <= sodium_emission_wl_nm <= filter_max_wl

    print("\n3. Comparing the tail's emission to the filter's range:")
    if not is_visible:
        print(f"   The sodium tail's light ({sodium_emission_wl_nm} nm) is outside the filter's range ({filter_min_wl}-{filter_max_wl} nm).")
        print("\n--- Conclusion ---")
        print("The filter will block the light from the sodium tail. Therefore, no tail would be visible.")
    else:
        # This case is not met, but included for completeness
        print(f"   The sodium tail's light ({sodium_emission_wl_nm} nm) is inside the filter's range ({filter_min_wl}-{filter_max_wl} nm).")
        print("\n--- Conclusion ---")
        print("The tail would be visible.")

if __name__ == '__main__':
    analyze_mercury_observation()