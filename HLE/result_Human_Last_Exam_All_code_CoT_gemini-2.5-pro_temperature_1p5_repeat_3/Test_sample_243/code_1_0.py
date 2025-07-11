import sys

def analyze_observation():
    """
    Analyzes the visibility of Mercury's sodium tail with a specific filter.
    """
    # Define observation parameters
    filter_center_wl = 486  # Filter center wavelength in nm
    filter_width = 1        # Filter width in nm

    # Define the emission lines of Mercury's sodium tail
    # The most prominent lines are the Sodium D-lines.
    sodium_d2_line = 589.0  # Wavelength in nm
    sodium_d1_line = 589.6  # Wavelength in nm

    # Calculate the filter's transmission range
    filter_lower_bound = filter_center_wl - (filter_width / 2.0)
    filter_upper_bound = filter_center_wl + (filter_width / 2.0)

    # Check if the sodium emission lines are within the filter's range
    is_d2_visible = filter_lower_bound <= sodium_d2_line <= filter_upper_bound
    is_d1_visible = filter_lower_bound <= sodium_d1_line <= filter_upper_bound

    # Print out the step-by-step analysis
    print("--- Analysis of Mercury's Tail Observation ---")
    print(f"Step 1: Define the filter's transmission range.")
    print(f"The filter is centered at {filter_center_wl} nm with a width of {filter_width} nm.")
    # The final equation is for the range calculation
    print(f"Equation for lower bound: {filter_center_wl} - ({filter_width} / 2.0) = {filter_lower_bound} nm")
    print(f"Equation for upper bound: {filter_center_wl} + ({filter_width} / 2.0) = {filter_upper_bound} nm")
    print(f"Resulting Filter Transmission Range: {filter_lower_bound} nm to {filter_upper_bound} nm\n")

    print(f"Step 2: Check if sodium emission lines fall within this range.")
    print(f"The sodium tail's primary color is yellow, from emission lines at ~{sodium_d2_line} nm and {sodium_d1_line} nm.")
    print(f"Is the {sodium_d2_line} nm line within the range [{filter_lower_bound}, {filter_upper_bound}]? {is_d2_visible}")
    print(f"Is the {sodium_d1_line} nm line within the range [{filter_lower_bound}, {filter_upper_bound}]? {is_d1_visible}\n")

    print("--- Conclusion ---")
    if not is_d2_visible and not is_d1_visible:
        print("The filter is designed to observe light at 486 nm (blue-green).")
        print(f"It completely blocks the yellow light from the sodium tail at ~{int(sodium_d2_line)} nm.")
        print("Therefore, no sodium tail would be visible with this specific filter.")
    else:
        print("The sodium tail would be visible through the filter.")

if __name__ == '__main__':
    analyze_observation()