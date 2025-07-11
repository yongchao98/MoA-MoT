def check_visibility():
    """
    Analyzes the visibility of Mercury's sodium tail based on emission
    and filter wavelengths.
    """
    # 1. Define the physical parameters from the problem.
    sodium_emission_wavelength = 589.0  # nm, for Sodium D-lines (yellow light)
    filter_center = 486.0               # nm, for the telescope's filter (blue light)
    filter_width = 1.0                  # nm, the bandwidth of the filter

    # 2. Calculate the filter's operational wavelength range.
    passband_min = filter_center - (filter_width / 2.0)
    passband_max = filter_center + (filter_width / 2.0)

    # 3. Check if the sodium emission is within the filter's range.
    # The final "equation" is this boolean check.
    is_sodium_light_visible = passband_min <= sodium_emission_wavelength <= passband_max

    # 4. Print the step-by-step reasoning and the conclusion.
    print("Analysis of Observation Feasibility")
    print("=" * 35)
    print(f"Wavelength of light from sodium tail: {sodium_emission_wavelength} nm")
    print(f"Wavelength passed by the telescope's filter: {filter_center} nm")
    print("\nThe filter allows a very narrow range of light to pass.")
    print(f"Filter Passband Calculation: {filter_center} nm +/- ({filter_width}/2) nm = [{passband_min}, {passband_max}] nm")
    print("\nVisibility Check:")
    print(f"Is the sodium light's wavelength ({sodium_emission_wavelength} nm) within the filter's range ({passband_min} nm to {passband_max} nm)?")
    print(f"Result: {is_sodium_light_visible}")
    print("=" * 35)

    if not is_sodium_light_visible:
        print("\nConclusion: The filter blocks the light from the sodium tail.")
        print("Because the sodium tail's light cannot pass through the filter, no tail would be visible.")
    else:
        # This case is not met by the problem's data.
        print("\nConclusion: The tail would be visible.")

check_visibility()