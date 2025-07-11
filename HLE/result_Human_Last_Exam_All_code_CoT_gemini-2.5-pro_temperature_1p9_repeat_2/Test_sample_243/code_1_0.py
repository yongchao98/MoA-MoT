import numpy as np

# Wavelength of the center of the bandpass filter in nanometers (nm)
filter_wavelength = 486
filter_color = "blue-green"

# Primary wavelength of emission for Mercury's sodium tail (nm)
# This is the average of the two sodium D-lines (589.0 and 589.6 nm)
sodium_tail_wavelength = 589.3
sodium_tail_color = "yellow"

# The filter is 1nm wide, so it passes light from 485.5 nm to 486.5 nm.
filter_passband_min = filter_wavelength - 0.5
filter_passband_max = filter_wavelength + 0.5

# Check if the sodium tail's light can pass through the filter
is_visible = filter_passband_min <= sodium_tail_wavelength <= filter_passband_max

print(f"Observation Attempt:")
print(f"Filter Center Wavelength: {filter_wavelength} nm (Color: {filter_color})")
print(f"Sodium Tail Emission Wavelength: {sodium_tail_wavelength} nm (Color: {sodium_tail_color})")
print("-" * 30)

if is_visible:
    print("Result: The sodium tail is visible through the filter.")
else:
    print(f"Result: The sodium tail is NOT visible.")
    print(f"Reason: The filter passes light from {filter_passband_min} nm to {filter_passband_max} nm.")
    print(f"The light from the sodium tail at {sodium_tail_wavelength} nm is outside this range and is blocked.")
    print("\nTherefore, no tail would be visible with this instrument setup.")
