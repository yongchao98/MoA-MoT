import numpy as np
# The 'colour-science' library is required for this script.
# You can install it by running: pip install colour-science
try:
    import colour
except ImportError:
    print("Please install the 'colour-science' library to run this script.")
    print("You can install it using pip: pip install colour-science")
    exit()

def solve():
    """
    Identifies which ROMM RGB colors are out of the sRGB gamut.
    """
    # Define the ROMM RGB colors to be tested, with their corresponding numbers.
    # ROMM RGB is technically 'ProPhoto RGB'.
    colors_to_test = {
        1: np.array([0, 0, 1]),
        2: np.array([0, 1, 0]),
        3: np.array([0, 0.5, 0.6]),
        4: np.array([0.4, 0.5, 0.6]),
        5: np.array([1, 1, 1])
    }

    # Define the input and output color spaces for the conversion.
    source_colorspace = colour.models.RGB_COLOURSPACE_PROPHOTO_RGB
    target_colorspace = colour.models.RGB_COLOURSPACE_sRGB

    out_of_gamut_indices = []

    print("Checking which ROMM RGB colors can be represented in sRGB...")
    print("-" * 60)

    # Iterate through each color, convert it, and check if it's in gamut.
    for number, romm_rgb in sorted(colors_to_test.items()):
        
        # Convert the color from ROMM RGB to sRGB
        srgb_rgb = colour.RGB_to_RGB(romm_rgb, source_colorspace, target_colorspace)

        # A color is in the sRGB gamut if all its components are between 0 and 1.
        is_in_gamut = np.all(srgb_rgb >= 0) and np.all(srgb_rgb <= 1)
        
        # Format the numbers for printing
        romm_r, romm_g, romm_b = romm_rgb
        srgb_r, srgb_g, srgb_b = srgb_rgb

        print(f"Color {number}: ROMM RGB({romm_r}, {romm_g}, {romm_b})")
        print(f"  -> Converted to sRGB: ({srgb_r:.4f}, {srgb_g:.4f}, {srgb_b:.4f})")
        
        if not is_in_gamut:
            result_text = "Cannot be represented in sRGB (Out of Gamut)."
            out_of_gamut_indices.append(number)
        else:
            result_text = "Can be represented in sRGB (In Gamut)."
            
        print(f"  -> Result: {result_text}\n")


    print("-" * 60)
    # Format the final answer as requested
    if not out_of_gamut_indices:
        final_answer = "none cannot"
    else:
        final_answer = ", ".join(map(str, sorted(out_of_gamut_indices)))
    
    print(f"Final Answer: The numbers of the colors that cannot be represented are:")
    print(final_answer)
    print(f"<<<{final_answer}>>>")

solve()