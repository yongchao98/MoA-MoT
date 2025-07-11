import sys

# This script checks which ROMM RGB colors are out of the sRGB gamut.
# It requires the 'colour-science' library for accurate color transformations.

try:
    import numpy as np
    import colour
except ImportError:
    print("This script requires the 'colour-science' and 'numpy' libraries.")
    print("Please install them by running: pip install colour-science numpy")
    sys.exit(1)

# A dictionary mapping the question number to its ROMM RGB color value
colors_to_test = {
    1: np.array([0, 0, 1]),      # Pure Blue
    2: np.array([0, 1, 0]),      # Pure Green
    3: np.array([0, 0.5, 0.6]),  # A cyan-like color
    4: np.array([0.4, 0.5, 0.6]),# A muted, grayish blue
    5: np.array([1, 1, 1]),      # White
}

# Define the source (ROMM RGB) and target (sRGB) color spaces using the library's models
romm_colorspace = colour.models.RGB_COLOURSPACE_PROPHOTO_RGB
srgb_colorspace = colour.models.RGB_COLOURSPACE_sRGB

print("Analyzing ROMM RGB to sRGB conversion for each color:")
print("-" * 70)

# Iterate through each color, convert it, and check if it's within the sRGB gamut
for number, romm_rgb_value in colors_to_test.items():
    # Convert the color from ROMM RGB to sRGB
    srgb_rgb_value = colour.RGB_to_RGB(romm_rgb_value, romm_colorspace, srgb_colorspace)

    # A color is out of gamut if any component is less than 0 or greater than 1
    is_out_of_gamut = np.any(srgb_rgb_value < 0) or np.any(srgb_rgb_value > 1)

    if is_out_of_gamut:
        status = "Out of sRGB Gamut"
    else:
        status = "In sRGB Gamut"

    # Show the input, output, and result of the check for each color.
    romm_str = f"({romm_rgb_value[0]:.1f}, {romm_rgb_value[1]:.1f}, {romm_rgb_value[2]:.1f})"
    srgb_str = f"({srgb_rgb_value[0]:.3f}, {srgb_rgb_value[1]:.3f}, {srgb_rgb_value[2]:.3f})"
    print(f"{number}) ROMM RGB {romm_str:<18} -> sRGB {srgb_str:<30} -> {status}")

print("-" * 70)