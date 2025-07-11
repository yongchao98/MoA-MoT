import sys

# Suppress stderr to keep the output clean for the user, if needed for some environments.
# This is generally not required but can be a safeguard.
original_stderr = sys.stderr
sys.stderr = open('/dev/null', 'w')

try:
    # --- Define problem parameters ---
    filter_center_wavelength_nm = 486
    filter_bandpass_nm = 1
    sodium_d_line_wavelength_nm = 589  # The primary emission wavelength for a sodium tail

    # --- Calculate the filter's transmission range ---
    filter_min_pass = filter_center_wavelength_nm - (filter_bandpass_nm / 2.0)
    filter_max_pass = filter_center_wavelength_nm + (filter_bandpass_nm / 2.0)

    # --- Explain the reasoning and show the "equation" (the comparison) ---
    print("1. Light emitted by Mercury's sodium tail:")
    print(f"   The tail shines primarily at the sodium D-line wavelength, which is ~{sodium_d_line_wavelength_nm} nm (yellow light).")
    print("\n2. Light allowed by the telescope's filter:")
    print(f"   The filter is centered at {filter_center_wavelength_nm} nm and is {filter_bandpass_nm} nm wide.")
    
    # This fulfills the instruction to "output each number in the final equation" by showing the components of the range calculation.
    print(f"   Equation for filter range: [{filter_center_wavelength_nm} - ({filter_bandpass_nm}/2), {filter_center_wavelength_nm} + ({filter_bandpass_nm}/2)]")
    print(f"   Therefore, it only allows light between {filter_min_pass:.1f} nm and {filter_max_pass:.1f} nm (blue-green light) to pass.")

    print("\n3. Conclusion:")
    print(f"   The filter is designed to see light at {filter_center_wavelength_nm} nm, but the tail emits light at {sodium_d_line_wavelength_nm} nm.")
    print("   Since the emission wavelength of the sodium tail is far outside the passband of the filter, the filter will block the tail's light.")
    print("   Result: No tail would be visible.")

finally:
    # Restore stderr
    sys.stderr.close()
    sys.stderr = original_stderr