import numpy as np
import colour

def check_romm_to_srgb_gamut():
    """
    Checks which ROMM RGB colors from a predefined list fall outside the sRGB gamut.
    """
    # Define the input ROMM RGB colors
    colors_to_check = {
        1: np.array([0, 0, 1.0]),     # Pure ROMM Blue
        2: np.array([0, 1.0, 0]),     # Pure ROMM Green
        3: np.array([0, 0.5, 0.6]),   # A cyan color
        4: np.array([0.4, 0.5, 0.6]), # A desaturated blue
        5: np.array([1.0, 1.0, 1.0])  # White point (D50)
    }

    # Define the source (ROMM RGB) and target (sRGB) color spaces
    romm_rgb_cs = colour.models.RGB_COLOURSPACE_PROPHOTO_RGB
    srgb_cs = colour.models.RGB_COLOURSPACE_sRGB

    out_of_gamut_numbers = []

    print("--- Checking ROMM RGB to sRGB Conversion ---")

    # Iterate through the colors, convert them, and check the gamut
    for number, romm_rgb in sorted(colors_to_check.items()):
        # Perform the color space conversion
        srgb_color = colour.RGB_to_RGB(romm_rgb, romm_rgb_cs, srgb_cs)

        # A color is out-of-gamut if any component is < 0 or > 1.
        # A small tolerance is used for floating-point inaccuracies.
        is_out_of_gamut = np.any((srgb_color < -1e-5) | (srgb_color > 1 + 1e-5))

        print(f"\nChecking Color #{number}: ROMM RGB({romm_rgb[0]}, {romm_rgb[1]}, {romm_rgb[2]})")
        print(f" -> Converted to sRGB: ({srgb_color[0]:.4f}, {srgb_color[1]:.4f}, {srgb_color[2]:.4f})")

        if is_out_of_gamut:
            out_of_gamut_numbers.append(number)
            print(" -> Result: CANNOT be represented in sRGB (Out of Gamut).")
        else:
            print(" -> Result: Can be represented in sRGB (In Gamut).")

    print("\n--- Final Answer ---")
    if not out_of_gamut_numbers:
        final_answer = "none cannot"
    else:
        final_answer = ", ".join(map(str, sorted(out_of_gamut_numbers)))

    print(f"The numbers of the colors that cannot be represented by an sRGB hex code are: {final_answer}")

if __name__ == '__main__':
    check_romm_to_srgb_gamut()
<<<1, 2, 3, 5>>>