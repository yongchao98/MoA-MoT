import numpy as np
try:
    import colour
except ImportError:
    print("Error: The 'colour-science' library is required.")
    print("Please install it using: pip install colour-science")
    exit()

def solve():
    """
    Determines which ROMM RGB colors are outside the sRGB gamut.
    """
    # Define the list of ROMM RGB colors to check, with their original numbers
    colors_romm = {
        1: np.array([0, 0, 1]),
        2: np.array([0, 1, 0]),
        3: np.array([0, 0.5, 0.6]),
        4: np.array([0.4, 0.5, 0.6]),
        5: np.array([1, 1, 1])
    }

    # Define the source (ROMM RGB) and target (sRGB) color spaces
    romm_colorspace = colour.RGB_COLOURSPACES['ProPhoto RGB']
    srgb_colorspace = colour.RGB_COLOURSPACES['sRGB']

    out_of_gamut_indices = []

    print("Checking which ROMM RGB values are outside the sRGB gamut...")
    print("-" * 60)

    # Iterate through each color, convert it, and check its gamut
    for index, romm_rgb in sorted(colors_romm.items()):
        # Convert from ROMM RGB to sRGB
        srgb_rgb = colour.RGB_to_RGB(romm_rgb, romm_colorspace, srgb_colorspace)

        # A color is representable in sRGB if all components are in [0, 1]
        is_in_gamut = np.all(srgb_rgb >= 0) and np.all(srgb_rgb <= 1)

        # Print the original values for the "equation"
        print(f"Case {index}: ROMM RGB({romm_rgb[0]}, {romm_rgb[1]}, {romm_rgb[2]})")
        print(f"  -> Converts to sRGB({srgb_rgb[0]:.4f}, {srgb_rgb[1]:.4f}, {srgb_rgb[2]:.4f})")

        if not is_in_gamut:
            print("  -> Result: Cannot be represented in sRGB (out of gamut).")
            out_of_gamut_indices.append(index)
        else:
            print("  -> Result: Can be represented in sRGB.")
        print("-" * 60)

    # Sort the indices in ascending order for the final output
    out_of_gamut_indices.sort()

    # Format the final answer string
    if not out_of_gamut_indices:
        final_answer_str = "none cannot"
    else:
        final_answer_str = ", ".join(map(str, out_of_gamut_indices))

    print("\nFinal Answer:")
    print("The numbers of the colors that cannot be represented are:")
    print(final_answer_str)

if __name__ == '__main__':
    solve()