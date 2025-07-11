import sys

def analyze_mercury_tail_visibility():
    """
    Analyzes the visibility of Mercury's sodium tail based on
    the emission wavelength and the observational filter used.
    """
    # Define the key wavelengths in nanometers (nm)
    # The dominant emission lines for sodium are the D-lines.
    sodium_emission_wavelength_d2 = 589.0
    sodium_emission_wavelength_d1 = 589.6
    
    # Define the properties of the observational filter
    filter_center_wavelength = 486.0
    filter_width = 1.0
    
    # Calculate the filter's passband (the range of wavelengths it allows through)
    filter_lower_bound = filter_center_wavelength - (filter_width / 2.0)
    filter_upper_bound = filter_center_wavelength + (filter_width / 2.0)
    
    # --- Print out the analysis ---
    print("Step 1: Define emission and observation parameters.")
    print(f"The primary emission wavelength of the sodium tail is ~{sodium_emission_wavelength_d2} nm (yellow light).")
    print(f"The telescope uses a filter centered at {filter_center_wavelength} nm (blue-green light) with a width of {filter_width} nm.")
    
    print("\nStep 2: Calculate the filter's allowed wavelength range.")
    # We use a simple equation: Range = Center ± (Width / 2)
    print(f"Filter Range Equation: {filter_center_wavelength} nm ± ({filter_width} nm / 2.0)")
    print(f"This results in a passband from {filter_lower_bound} nm to {filter_upper_bound} nm.")
    
    # Check if the sodium emission falls within the filter's range
    is_visible = (filter_lower_bound <= sodium_emission_wavelength_d2 <= filter_upper_bound) or \
                 (filter_lower_bound <= sodium_emission_wavelength_d1 <= filter_upper_bound)
                 
    print("\nStep 3: Compare emission wavelength to the filter's range.")
    print(f"Is the sodium emission wavelength ({sodium_emission_wavelength_d2} nm) within the filter's range ({filter_lower_bound}-{filter_upper_bound} nm)? {is_visible}")
    
    print("\n--- Conclusion ---")
    if is_visible:
        print("The sodium tail would be visible.")
    else:
        print("The filter blocks the light from the sodium tail.")
        print("Therefore, no tail would be visible.")

# Run the analysis
analyze_mercury_tail_visibility()
