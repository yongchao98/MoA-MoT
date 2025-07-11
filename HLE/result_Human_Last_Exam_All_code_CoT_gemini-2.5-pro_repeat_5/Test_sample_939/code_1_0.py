import numpy as np

def solve():
    """
    Determines which ROMM RGB color values cannot be represented by an sRGB hex code.
    """

    # Standard matrix to convert linear ProPhoto RGB (D50) to CIE XYZ (D50).
    M_ROMM_TO_XYZ = np.array([
        [0.7976749, 0.1351917, 0.0313534],
        [0.2880402, 0.7118741, 0.0000857],
        [0.0,       0.0,       0.8252100]
    ])

    # Standard Bradford chromatic adaptation matrix to convert from D50 illuminant to D65.
    M_CAT_D50_TO_D65 = np.array([
        [0.9555766, -0.0230393, 0.0631636],
        [-0.0282895, 1.0099416, 0.0210077],
        [0.0122982, -0.0204830, 1.3299098]
    ])

    # Standard matrix to convert CIE XYZ (D65) to linear sRGB.
    M_XYZ_TO_sRGB = np.array([
        [ 3.24096994, -1.53738335, -0.49861076],
        [-0.96924363,  1.8759675,   0.04155506],
        [ 0.05563001, -0.20397696,  1.05697151]
    ])

    # Combine matrices for a single-step conversion from linear ROMM RGB to linear sRGB.
    M_conversion = M_XYZ_TO_sRGB @ M_CAT_D50_TO_D65 @ M_ROMM_TO_XYZ

    colors = [
        np.array([0, 0, 1]),
        np.array([0, 1, 0]),
        np.array([0, 0.5, 0.6]),
        np.array([0.4, 0.5, 0.6]),
        np.array([1, 1, 1])
    ]

    out_of_gamut_indices = []
    
    # Define the gamma for ROMM RGB
    gamma = 1.8

    for i, color in enumerate(colors):
        color_num = i + 1
        print(f"Checking color {color_num}: RGB({color[0]}, {color[1]}, {color[2]})")

        # 1. Linearize the ROMM RGB values
        linear_romm = np.power(color, gamma)
        print(f"  - Linear ROMM RGB: ({linear_romm[0]:.4f}, {linear_romm[1]:.4f}, {linear_romm[2]:.4f})")
        
        # 2. Convert to linear sRGB
        linear_srgb = M_conversion @ linear_romm
        
        r, g, b = linear_srgb[0], linear_srgb[1], linear_srgb[2]
        
        print(f"  - Converted to Linear sRGB: R={r:.4f}, G={g:.4f}, B={b:.4f}")

        # 3. Check if the color is in the sRGB gamut
        # A small tolerance is used for floating point comparisons.
        tolerance = 1e-9
        if r < -tolerance or r > 1 + tolerance or \
           g < -tolerance or g > 1 + tolerance or \
           b < -tolerance or b > 1 + tolerance:
            print("  - Result: Cannot be represented in sRGB (out of gamut).")
            out_of_gamut_indices.append(str(color_num))
        else:
            print("  - Result: Can be represented in sRGB (in gamut).")
        print("-" * 20)

    if not out_of_gamut_indices:
        answer = "none cannot"
    else:
        answer = ", ".join(out_of_gamut_indices)

    print(f"\nThe numbers of the colors that cannot be represented are: {answer}")
    print(f"<<<{answer}>>>")

solve()