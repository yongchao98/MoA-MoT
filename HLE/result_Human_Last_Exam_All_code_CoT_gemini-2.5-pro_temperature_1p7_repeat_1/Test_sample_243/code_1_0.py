import sys

def analyze_observation():
    """
    Analyzes the visibility of Mercury's sodium tail with a specific filter.
    """
    # Wavelength of Sodium D-line emission (primary emission from the sodium tail)
    sodium_emission_nm = 589.0

    # Properties of the bandpass filter from the problem
    filter_center_nm = 486.0
    filter_width_nm = 1.0

    # Calculate the filter's passband range
    filter_min_nm = filter_center_nm - (filter_width_nm / 2.0)
    filter_max_nm = filter_center_nm + (filter_width_nm / 2.0)

    # Explain the setup
    print(f"Step 1: Identify the wavelength of light from Mercury's sodium tail.")
    print(f"Sodium atoms in the tail primarily emit yellow light at ~{sodium_emission_nm} nm.\n")

    print(f"Step 2: Identify the properties of the observation filter.")
    print(f"The filter is centered at {filter_center_nm} nm (blue light).")
    print(f"The filter width is {filter_width_nm} nm.\n")

    print(f"Step 3: Calculate the range of light the filter allows to pass.")
    # This is the "equation" part requested by the prompt
    print(f"The filter passband is calculated as: {filter_center_nm} +/- ({filter_width_nm} / 2)")
    print(f"This means it only allows light from {filter_min_nm} nm to {filter_max_nm} nm to pass through.\n")

    print("Step 4: Compare the sodium tail's light to the filter's range.")
    print(f"Is the sodium emission ({sodium_emission_nm} nm) within the filter's range ({filter_min_nm} nm to {filter_max_nm} nm)?")

    # Check if the sodium emission is visible through the filter
    if filter_min_nm <= sodium_emission_nm <= filter_max_nm:
        print("\nConclusion: Yes, the tail would be visible.")
    else:
        print("\nConclusion: No. The light from the sodium tail is at a different wavelength than what the filter allows.")
        print("The filter will block the light from the sodium tail.")
        print("Therefore, no tail would be visible.")

if __name__ == '__main__':
    analyze_observation()