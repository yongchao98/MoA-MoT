import math

def analyze_observation():
    """
    Analyzes the visibility of Mercury's sodium tail based on a specific filter.
    """

    # --- Define Known Physical Constants ---
    # Wavelength of Mercury's sodium tail emission in nanometers (nm)
    sodium_emission_wavelength = 589.0
    # Color corresponding to the sodium emission
    sodium_emission_color = "yellow"

    # Wavelength of the filter center in nanometers (nm)
    filter_center_wavelength = 486.0
    # Bandwidth of the filter in nanometers (nm)
    filter_bandwidth = 1.0

    # --- Step 1: Explain the physics of the tail length ---
    print("--- Analysis of Mercury's Tail ---")
    print("Mercury's tail length depends on its distance from the Sun:")
    print(" - At Perihelion (closest), solar radiation is strongest, creating a LONGER tail.")
    print(" - At Aphelion (farthest), solar radiation is weakest, creating a SHORTER tail.")
    print("-" * 20)

    # --- Step 2: Determine what the filter lets through ---
    filter_lower_bound = filter_center_wavelength - (filter_bandwidth / 2.0)
    filter_upper_bound = filter_center_wavelength + (filter_bandwidth / 2.0)

    print("--- Analysis of the Observation ---")
    print(f"The instrument is using a filter that only allows light from {filter_lower_bound} nm to {filter_upper_bound} nm to pass through.")
    print(f"The sodium tail emits light primarily at {sodium_emission_wavelength} nm, which has a {sodium_emission_color} color.")
    print("-" * 20)

    # --- Step 3: Check if the tail's light passes through the filter ---
    print("--- Conclusion ---")
    # Using 'is not' is more Pythonic for checking against a range
    is_visible = filter_lower_bound <= sodium_emission_wavelength <= filter_upper_bound

    if is_visible:
        print(f"The sodium tail's light at {sodium_emission_wavelength} nm would be visible through the filter.")
    else:
        print(f"The filter (passband {filter_lower_bound}-{filter_upper_bound} nm) BLOCKS the light from the sodium tail (emission at {sodium_emission_wavelength} nm).")
        print("Therefore, no tail would be visible.")


analyze_observation()