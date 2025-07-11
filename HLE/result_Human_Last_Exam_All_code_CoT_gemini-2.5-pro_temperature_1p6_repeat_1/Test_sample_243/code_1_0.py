import sys

def analyze_mercury_tail_observation():
    """
    Analyzes the visibility of Mercury's sodium tail given specific
    observational parameters.
    """
    # --- Step 1: Define parameters from the problem ---
    # The filter is centered at 486 nm
    filter_center_wavelength = 486
    # The filter has a bandwidth of 1 nm
    filter_bandwidth = 1
    # Mercury's tail is primarily sodium, which emits strongly at the sodium D-lines
    sodium_emission_wavelength = 589

    print("Analyzing the visibility of Mercury's sodium tail...")
    print("-" * 50)
    print(f"Observation Filter Wavelength: {filter_center_wavelength} nm")
    print(f"Observation Filter Bandwidth: {filter_bandwidth} nm")
    print(f"Source (Sodium Tail) Emission Wavelength: {sodium_emission_wavelength} nm")
    print("-" * 50)

    # --- Step 2: Calculate the filter's transmission window ---
    # The equation for the filter's range is [center - (bandwidth/2), center + (bandwidth/2)]
    lower_bound = filter_center_wavelength - (filter_bandwidth / 2)
    upper_bound = filter_center_wavelength + (filter_bandwidth / 2)

    # --- Step 3: Output the check being performed ---
    # Fulfills the requirement to show each number in the "equation" or logical check.
    print("Checking if the sodium emission falls within the filter's range.")
    print(f"Filter transmission equation: Light is passed if wavelength 'w' satisfies:")
    print(f"{lower_bound} <= w <= {upper_bound}")
    print(f"Is {sodium_emission_wavelength} (sodium) between {lower_bound} and {upper_bound}?")

    # --- Step 4: Determine visibility and provide a final conclusion ---
    is_visible = lower_bound <= sodium_emission_wavelength <= upper_bound

    if is_visible:
        # This path is logically possible but incorrect for the given numbers.
        # It is included for a complete logical flow.
        print("\nConclusion: The sodium emission IS within the filter's range.")
        print("The tail would be visible.")
    else:
        print("\nConclusion: The sodium emission IS NOT within the filter's range.")
        print("The 589 nm light from the sodium tail is blocked by the 486 nm filter.")
        print("Therefore, no tail would be visible.")

if __name__ == '__main__':
    analyze_mercury_tail_observation()
