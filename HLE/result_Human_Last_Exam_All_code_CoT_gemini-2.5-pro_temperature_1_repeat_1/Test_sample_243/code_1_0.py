def analyze_mercury_observation(filter_center_nm, filter_width_nm, tail_emission_nm):
    """
    Analyzes the visibility of Mercury's sodium tail based on filter properties.
    """
    # Calculate the range of wavelengths the filter allows to pass
    filter_min_nm = filter_center_nm - (filter_width_nm / 2)
    filter_max_nm = filter_center_nm + (filter_width_nm / 2)

    print("--- Analysis of Mercury's Sodium Tail Observation ---")
    print(f"Step 1: Analyzing the effect of the orbital position.")
    print("The tail's length is dependent on solar radiation. It is longer at perihelion (closer to Sun) and shorter at aphelion (farther from Sun).\n")

    print(f"Step 2: Analyzing the effect of the observation filter.")
    print(f"The observation uses a filter centered at {filter_center_nm} nm with a width of {filter_width_nm} nm.")
    print(f"This means the filter only passes light between {filter_min_nm:.1f} nm and {filter_max_nm:.1f} nm.\n")

    print(f"The sodium tail's characteristic emission is at approximately {tail_emission_nm} nm (yellow light).")

    # Check if the tail's emission is within the filter's range
    if filter_min_nm <= tail_emission_nm <= filter_max_nm:
        print("\nConclusion: The tail's light is within the filter's range and would be visible.")
    else:
        print(f"\nConclusion: The tail's {tail_emission_nm} nm light is outside the filter's range and will be blocked.")
        print("Therefore, no tail would be visible through this instrument.")

# Parameters from the problem
filter_center = 486
filter_width = 1
sodium_emission = 589

# Run the analysis
analyze_mercury_observation(filter_center, filter_width, sodium_emission)