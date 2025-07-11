# Define the physical properties from the problem description
sodium_emission_wavelength = 589.0  # Dominant wavelength for sodium in nm (yellow)
filter_center_wavelength = 486.0    # Center of the bandpass filter in nm (blue-green)
filter_width = 1.0                  # Width of the bandpass filter in nm

# Calculate the wavelength range the filter allows to pass through
filter_min = filter_center_wavelength - (filter_width / 2.0)
filter_max = filter_center_wavelength + (filter_width / 2.0)

# Print the parameters for clarity
print(f"The observation is being made with a filter centered at {filter_center_wavelength} nm.")
print(f"This filter only allows light between {filter_min:.1f} nm and {filter_max:.1f} nm to pass through.")
print(f"Mercury's sodium tail primarily emits light at the sodium D-lines, around {sodium_emission_wavelength} nm.")

# Check if the sodium emission wavelength is within the filter's range
if filter_min <= sodium_emission_wavelength <= filter_max:
    print("\nConclusion: The sodium tail's light is within the filter's range and would be visible.")
else:
    print(f"\nConclusion: The sodium tail's light at {sodium_emission_wavelength} nm is outside the filter's range of {filter_min:.1f}-{filter_max:.1f} nm.")
    print("Therefore, the filter blocks the light from the sodium tail, and no tail would be visible.")
