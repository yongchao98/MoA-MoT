import sys

def analyze_mercury_tail_visibility():
    """
    Analyzes the visibility of Mercury's sodium tail based on
    the properties of an observational filter.
    """

    # --- Step 1: Define properties of the observational filter ---
    # The filter is centered at 486 nm with a 1 nm bandpass.
    filter_center_nm = 486
    filter_bandwidth_nm = 1
    filter_min_pass_nm = filter_center_nm - (filter_bandwidth_nm / 2.0)
    filter_max_pass_nm = filter_center_nm + (filter_bandwidth_nm / 2.0)

    # --- Step 2: Define properties of Mercury's sodium tail ---
    # The tail's glow comes from sodium atoms emitting light at specific wavelengths.
    # The most prominent emission lines for neutral sodium (Na) are the D-lines, which are yellow.
    sodium_d_line_1_nm = 589.0
    sodium_d_line_2_nm = 589.6

    # --- Step 3: Print the analysis ---
    print("Analysis of Mercury's Tail Observation")
    print("=" * 40)
    print(f"The telescope is using a filter that allows light from {filter_min_pass_nm} nm to {filter_max_pass_nm} nm.")
    print(f"The primary visible light from the sodium tail is at the sodium D-lines: {sodium_d_line_1_nm} nm and {sodium_d_line_2_nm} nm.")
    print("\nComparing the filter's range to the sodium emission:")

    # --- Step 4: Draw a conclusion based on the comparison ---
    # Check if the sodium light is within the filter's range.
    is_visible_d1 = filter_min_pass_nm <= sodium_d_line_1_nm <= filter_max_pass_nm
    is_visible_d2 = filter_min_pass_nm <= sodium_d_line_2_nm <= filter_max_pass_nm

    if not is_visible_d1 and not is_visible_d2:
        print(f"The filter (centered at {filter_center_nm} nm) does NOT allow light from the sodium lines (~589 nm) to pass.")
        print("\nFinal Conclusion: Because the filter blocks the characteristic light of the sodium tail, no tail would be visible.")
        # Note: The filter at 486 nm is actually for observing hydrogen (the H-beta line is at 486.1 nm), not sodium.
    else:
        # This case is not physically correct given the parameters.
        print("The filter allows sodium light to pass.")
        print("\nFinal Conclusion: The tail would be visible and its length would vary.")

# Execute the analysis
analyze_mercury_tail_visibility()
