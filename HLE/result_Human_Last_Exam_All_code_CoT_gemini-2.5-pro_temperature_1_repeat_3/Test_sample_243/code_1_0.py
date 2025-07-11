import sys

# Define the wavelengths and filter properties
sodium_d_line_wavelength = 589.0  # nm, the primary emission wavelength for sodium
filter_center_wavelength = 486.0  # nm
filter_width = 1.0  # nm

# Calculate the filter's passband range
filter_min = filter_center_wavelength - filter_width / 2.0
filter_max = filter_center_wavelength + filter_width / 2.0

# Print the initial parameters
print(f"Mercury's sodium tail emits light primarily at: {sodium_d_line_wavelength} nm (Yellow)")
print(f"The telescope's filter is centered at: {filter_center_wavelength} nm")
print(f"The filter allows light in the range: {filter_min} nm to {filter_max} nm")
print("-" * 30)

# Check if the sodium emission is visible through the filter
if filter_min <= sodium_d_line_wavelength <= filter_max:
    print("Result: The sodium emission is WITHIN the filter's range.")
    print("Conclusion: The tail would be visible.")
else:
    print("Result: The sodium emission is OUTSIDE the filter's range.")
    print("Conclusion: The filter blocks the light from the sodium tail. Therefore, no tail would be visible.")
