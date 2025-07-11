# Step 1: Define the emission wavelengths for the sodium tail.
# The most prominent emission lines for sodium are the D-lines.
sodium_d1_wavelength_nm = 589.6
sodium_d2_wavelength_nm = 589.0

# Step 2: Define the properties of the optical filter.
filter_center_nm = 486
filter_width_nm = 1
filter_min_nm = filter_center_nm - (filter_width_nm / 2)
filter_max_nm = filter_center_nm + (filter_width_nm / 2)

# Step 3: Check if the sodium emission lines fall within the filter's passband.
is_d1_visible = filter_min_nm <= sodium_d1_wavelength_nm <= filter_max_nm
is_d2_visible = filter_min_nm <= sodium_d2_wavelength_nm <= filter_max_nm

# Step 4: Print the analysis and conclusion.
print("Analysis of Mercury's Sodium Tail Visibility")
print("---------------------------------------------")
print(f"The telescope's filter is centered at {filter_center_nm}nm and passes light from {filter_min_nm}nm to {filter_max_nm}nm.")
print(f"The sodium tail primarily emits yellow light at approximately {sodium_d2_wavelength_nm}nm and {sodium_d1_wavelength_nm}nm.")
print(f"Checking if the {sodium_d2_wavelength_nm}nm sodium line is visible through the filter: {is_d2_visible}")
print(f"Checking if the {sodium_d1_wavelength_nm}nm sodium line is visible through the filter: {is_d1_visible}")
print("\nConclusion:")
if not is_d1_visible and not is_d2_visible:
    print("The filter completely blocks the light from the sodium tail.")
    print("Therefore, no tail would be visible, regardless of Mercury's position (aphelion or perihelion).")
else:
    print("The sodium tail would be visible through the filter.")
