def check_tail_visibility():
    """
    Analyzes the visibility of Mercury's sodium tail based on emission lines
    and a specified optical filter.
    """
    # Key physical constants in nanometers (nm)
    # The sodium tail is primarily visible due to the Sodium D-lines.
    sodium_d_line_wavelength_nm = 589.0  # We can use the stronger D2 line for this example.

    # Properties of the observational filter in nm
    filter_center_nm = 486
    filter_width_nm = 1

    # Calculate the passband of the filter
    # The filter allows light to pass from half the width below the center
    # to half the width above the center.
    filter_min_pass_nm = filter_center_nm - (filter_width_nm / 2)
    filter_max_pass_nm = filter_center_nm + (filter_width_nm / 2)

    # Output the problem setup and the logical steps
    print("--- Analyzing Visibility of Mercury's Sodium Tail ---")
    print(f"The tail's primary emission wavelength (Sodium D-line) is: {sodium_d_line_wavelength_nm} nm")
    print(f"The camera's filter is centered at: {filter_center_nm} nm")
    print(f"The camera's filter has a width of: {filter_width_nm} nm")
    print("-" * 20)

    # Perform and print the visibility check, showing all numbers involved
    print("To be visible, the tail's emission must be within the filter's passband.")
    print(f"Filter Passband Calculation: {filter_center_nm} ± ({filter_width_nm} / 2) = {filter_min_pass_nm} nm to {filter_max_pass_nm} nm")
    print("-" * 20)

    # Final check and conclusion
    is_visible = filter_min_pass_nm <= sodium_d_line_wavelength_nm <= filter_max_pass_nm
    
    print(f"Is {sodium_d_line_wavelength_nm} nm within the range [{filter_min_pass_nm}, {filter_max_pass_nm}]? -> {is_visible}")

    if is_visible:
        print("\nConclusion: The sodium tail WOULD be visible through the filter.")
    else:
        print("\nConclusion: The filter BLOCKS the light from the sodium tail.")
        print("Therefore, no tail would be visible, regardless of its length or brightness.")

check_tail_visibility()