import math

# Step 1: Define the characteristic emission wavelength of the sodium tail.
# The bright yellow light from sodium (like in streetlights) is from its D-lines.
sodium_emission_wavelength_nm = 589.0

# Step 2: Define the parameters of the observation equipment.
filter_center_nm = 486
filter_width_nm = 1

# Step 3: Calculate the filter's passband range.
half_width = filter_width_nm / 2
filter_min_nm = filter_center_nm - half_width
filter_max_nm = filter_center_nm + half_width

# Step 4: Compare the sodium emission to the filter's passband.
is_visible = filter_min_nm <= sodium_emission_wavelength_nm <= filter_max_nm

# Step 5: Print the analysis and conclusion.
print("Analysis of Mercury's Sodium Tail Visibility")
print("="*45)
print(f"The sodium that makes up Mercury's tail emits light at ~{sodium_emission_wavelength_nm} nm.")
print(f"The telescope's filter is centered at {filter_center_nm} nm with a {filter_width_nm} nm width.")
print(f"The filter's passband equation is: filter_center ± (filter_width / 2)")
print(f"The calculation is: {filter_center_nm} ± ({filter_width_nm} / 2), which gives a range of {filter_min_nm} nm to {filter_max_nm} nm.")
print(f"\nChecking if sodium's light ({sodium_emission_wavelength_nm} nm) falls within the filter's range [{filter_min_nm}, {filter_max_nm}] nm...")

if is_visible:
    print("\nResult: The sodium emission is WITHIN the filter's range.")
    # This case is physically incorrect but included for logical completeness.
else:
    print("\nResult: The sodium emission is OUTSIDE the filter's range.")
    print("Conclusion: Because the filter blocks the light from the sodium tail, no tail would be visible.")
    print("The differences between aphelion and perihelion are irrelevant for this observation.")
