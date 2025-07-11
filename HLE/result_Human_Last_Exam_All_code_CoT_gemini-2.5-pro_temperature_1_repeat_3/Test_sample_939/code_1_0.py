import numpy as np

# The user might need to install the 'colour-science' library first.
# This can be done by running: pip install colour-science
try:
    import colour
except ImportError:
    print("Error: The 'colour-science' library is required.")
    print("Please install it by running: pip install colour-science")
    exit()

def check_romm_to_srgb_gamut():
    """
    Checks which ROMM RGB colors are outside the sRGB gamut.
    """
    # Define the colors to be tested, with their corresponding numbers
    colors_to_test = {
        1: np.array([0., 0., 1.]),   # Pure ROMM Blue
        2: np.array([0., 1., 0.]),   # Pure ROMM Green
        3: np.array([0., 0.5, 0.6]), # A saturated Cyan
        4: np.array([0.4, 0.5, 0.6]),# A less saturated color
        5: np.array([1., 1., 1.])    # White point
    }

    out_of_gamut_numbers = []

    # Define the input and output color spaces using the colour-science library
    romm_colorspace = colour.models.PROPHOTO_RGB_COLOURSPACE
    srgb_colorspace = colour.models.sRGB_COLOURSPACE

    print("Converting ROMM RGB values to sRGB to check for out-of-gamut colors...")
    print("-" * 70)
    print(f"{'ID':<5}{'Input ROMM RGB':<20}{'Resulting sRGB':<30}{'In Gamut?'}")
    print("-" * 70)

    for number, romm_rgb in colors_to_test.items():
        # Perform the conversion
        srgb_color = colour.RGB_to_RGB(romm_rgb, romm_colorspace, srgb_colorspace)

        # A color is out-of-gamut if any component is < 0 or > 1.
        # We use a small tolerance for floating point math.
        is_in_gamut = np.all((srgb_color >= -1e-7) & (srgb_color <= 1.0 + 1e-7))

        if not is_in_gamut:
            out_of_gamut_numbers.append(str(number))

        # Print the results for each color
        romm_str = f"({romm_rgb[0]:.1f}, {romm_rgb[1]:.1f}, {romm_rgb[2]:.1f})"
        
        # In the final equation, we show the converted sRGB values
        srgb_str = f"({srgb_color[0]:.4f}, {srgb_color[1]:.4f}, {srgb_color[2]:.4f})"
        
        print(f"{number:<5}{romm_str:<20}{srgb_str:<30}{is_in_gamut}")

    print("-" * 70)
    
    # Phrase the final answer as requested
    if not out_of_gamut_numbers:
        final_answer = "none cannot"
    else:
        final_answer = ", ".join(sorted(out_of_gamut_numbers))

    print(f"\nThe numbers of the colors that cannot be represented are:\n{final_answer}")
    
    # Return the final answer in the specified format
    print(f"\n<<< {final_answer} >>>")


if __name__ == '__main__':
    check_romm_to_srgb_gamut()