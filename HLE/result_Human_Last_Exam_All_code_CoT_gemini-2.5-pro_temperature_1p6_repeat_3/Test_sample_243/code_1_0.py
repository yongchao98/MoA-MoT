import math

# Define the filter parameters
filter_center_wl = 486  # in nanometers
filter_width = 1.0      # in nanometers

# Calculate the range of the filter
filter_min_wl = filter_center_wl - (filter_width / 2.0)
filter_max_wl = filter_center_wl + (filter_width / 2.0)

# Define the prominent emission wavelengths of Sodium (D-lines)
sodium_d1_wl = 589.6 # in nanometers
sodium_d2_wl = 589.0 # in nanometers

# Explain the analysis
print(f"Step 1: Analyzing the observational setup.")
print(f"The telescope uses a filter centered at {filter_center_wl} nm with a {filter_width} nm width.")
print(f"This filter only allows light between {filter_min_wl} nm and {filter_max_wl} nm to pass through.")
print(f"The wavelength of {filter_center_wl} nm corresponds to a blue-green color in the visible spectrum.\n")

print(f"Step 2: Analyzing the light source (Mercury's tail).")
print(f"Mercury's tail is composed of sodium atoms, which primarily emit yellow-orange light at two wavelengths: {sodium_d2_wl} nm and {sodium_d1_wl} nm.\n")

print(f"Step 3: Comparing the light source to the filter.")
# Check if the sodium light can pass through the filter
can_see_d1 = filter_min_wl <= sodium_d1_wl <= filter_max_wl
can_see_d2 = filter_min_wl <= sodium_d2_wl <= filter_max_wl

print(f"Checking if Sodium D1 line ({sodium_d1_wl} nm) is within the filter range ({filter_min_wl} nm - {filter_max_wl} nm): {can_see_d1}")
print(f"Checking if Sodium D2 line ({sodium_d2_wl} nm) is within the filter range ({filter_min_wl} nm - {filter_max_wl} nm): {can_see_d2}\n")


print("Step 4: Conclusion.")
if not can_see_d1 and not can_see_d2:
    print("The filter completely blocks the yellow light from the sodium tail.")
    print("Therefore, no tail would be visible regardless of Mercury's position (aphelion or perihelion).")
    print("The correct answer choice describes this outcome.")
else:
    print("The sodium tail would be visible.")
