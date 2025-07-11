import numpy as np
try:
    from skimage.color import convert_colorspace
except ImportError:
    print("Please install the required libraries: pip install scikit-image numpy")
    exit()

def check_romm_to_srgb_gamut():
    """
    Converts a list of ROMM RGB colors to sRGB and checks if they are within the sRGB gamut.
    """
    # Define the ROMM RGB colors from the problem description
    colors = {
        1: np.array([0, 0, 1]),    # Pure ROMM Blue
        2: np.array([0, 1, 0]),    # Pure ROMM Green
        3: np.array([0, 0.5, 0.6]),
        4: np.array([0.4, 0.5, 0.6]),
        5: np.array([1, 1, 1]),    # White
    }

    out_of_gamut_colors = []

    print("--- Checking ROMM RGB to sRGB Conversion ---")
    
    for number, romm_rgb in colors.items():
        # The conversion function expects an image-like array, so we reshape the color
        # from a simple array [R, G, B] to a (1, 1, 3) array.
        romm_rgb_reshaped = romm_rgb.reshape((1, 1, 3))

        # Perform the color space conversion from ProPhoto (ROMM) to sRGB
        srgb_reshaped = convert_colorspace(romm_rgb_reshaped, 'ProPhoto RGB', 'sRGB')

        # Squeeze the result back to a simple [R, G, B] array
        srgb = srgb_reshaped.squeeze()
        r, g, b = srgb[0], srgb[1], srgb[2]

        print(f"\nColor {number}: ROMM RGB({romm_rgb[0]}, {romm_rgb[1]}, {romm_rgb[2]})")
        print(f"  -> Converted to sRGB({r:.4f}, {g:.4f}, {b:.4f})")

        # Check if any sRGB component is outside the valid [0, 1] range.
        is_in_gamut = (0 <= r <= 1) and (0 <= g <= 1) and (0 <= b <= 1)

        if not is_in_gamut:
            status = "Out of sRGB gamut. CANNOT be represented."
            out_of_gamut_colors.append(number)
        else:
            status = "In sRGB gamut. Can be represented."
            
        print(f"  -> Status: {status}")

    print("\n----------------------------------------------")
    
    # Format the final answer
    if not out_of_gamut_colors:
        final_answer = "none cannot"
    else:
        out_of_gamut_colors.sort()
        final_answer = ", ".join(map(str, out_of_gamut_colors))

    print(f"The numbers of the colors that cannot be represented are: {final_answer}")
    print(f"<<<{final_answer}>>>")


if __name__ == '__main__':
    check_romm_to_srgb_gamut()