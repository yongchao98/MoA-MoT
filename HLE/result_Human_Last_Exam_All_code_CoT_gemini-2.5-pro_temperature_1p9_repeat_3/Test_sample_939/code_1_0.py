import numpy as np

def solve():
    """
    Determines which ROMM RGB colors are outside the sRGB gamut.
    """
    # Matrix to convert from linear ROMM RGB (D50) to linear sRGB (D65)
    # This matrix accounts for both primary rotation and chromatic adaptation.
    PROPHOTO_TO_SRGB_MATRIX = np.array([
        [ 2.0413690, -0.5649464, -0.3446944],
        [-0.9692660,  1.8760108,  0.0415560],
        [ 0.0134474, -0.1183897,  1.0154096]
    ])

    # List of ROMM RGB colors to test
    colors = {
        1: np.array([0, 0, 1]),
        2: np.array([0, 1, 0]),
        3: np.array([0, 0.5, 0.6]),
        4: np.array([0.4, 0.5, 0.6]),
        5: np.array([1, 1, 1])
    }

    out_of_gamut_indices = []

    print("Checking which ROMM RGB colors can be represented in sRGB...\n")

    for i in sorted(colors.keys()):
        romm_rgb = colors[i]
        
        # Convert ROMM RGB to linear sRGB
        srgb_linear = PROPHOTO_TO_SRGB_MATRIX @ romm_rgb
        
        # Check if the converted sRGB values are in the gamut [0, 1]
        is_in_gamut = np.all(srgb_linear >= 0) and np.all(srgb_linear <= 1)
        
        print(f"Color {i}: ROMM RGB({romm_rgb[0]}, {romm_rgb[1]}, {romm_rgb[2]})")
        print(f" -> Converted to sRGB({srgb_linear[0]:.4f}, {srgb_linear[1]:.4f}, {srgb_linear[2]:.4f})")

        if not is_in_gamut:
            out_of_gamut_indices.append(i)
            print(" -> Result: Cannot be represented (Out of sRGB gamut)\n")
        else:
            print(" -> Result: Can be represented (In sRGB gamut)\n")

    if not out_of_gamut_indices:
        answer = "none cannot"
    else:
        answer = ", ".join(map(str, sorted(out_of_gamut_indices)))
        
    print(f"Final Answer: The numbers of the colors that cannot be represented are:")
    print(answer)
    print(f"<<<{answer}>>>")

solve()