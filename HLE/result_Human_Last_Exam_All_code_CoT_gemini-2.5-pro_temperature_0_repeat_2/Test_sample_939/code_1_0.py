import numpy as np

def solve():
    """
    Determines which ROMM RGB colors are outside the sRGB gamut.
    """

    # --- 1. Define Conversion Matrices ---

    # Matrix to convert from ProPhoto (ROMM) RGB to CIE XYZ with a D50 white point.
    M_ROMM_to_XYZ_D50 = np.array([
        [0.7976749, 0.1351917, 0.0313534],
        [0.2880402, 0.7118741, 0.0000857],
        [0.0000000, 0.0000000, 0.8252100]
    ])

    # Bradford chromatic adaptation matrix to convert from D50 white point to D65.
    M_D50_to_D65 = np.array([
        [ 0.9555766, -0.0230393,  0.0631636],
        [-0.0282895,  1.0099416,  0.0210077],
        [ 0.0122982, -0.0204830,  1.3299098]
    ])

    # Matrix to convert from CIE XYZ (D65) to linear sRGB.
    M_XYZ_D65_to_sRGB = np.array([
        [ 3.2404542, -1.5371385, -0.4985314],
        [-0.9692660,  1.8760108,  0.0415560],
        [ 0.0556434, -0.2040259,  1.0572252]
    ])

    # Combine matrices for a single transformation from linear ROMM (D50) to linear sRGB (D65).
    M_final_transform = M_XYZ_D65_to_sRGB @ M_D50_to_D65 @ M_ROMM_to_XYZ_D50

    # --- 2. Define Helper Function and Input Colors ---

    def delinearize_romm_component(c):
        """Applies the inverse ROMM RGB transfer function to a single component."""
        # ProPhoto RGB (ROMM RGB) transfer function EOTF^-1
        # V_linear = C' ^ 1.8 for C' >= 16/512
        # V_linear = C' / 16 for C' < 16/512
        E = 16.0 / 512.0  # 0.03125
        if c >= E:
            return c ** 1.8
        else:
            return c / 16.0

    colors = [
        (0, 0, 1),
        (0, 1, 0),
        (0, 0.5, 0.6),
        (0.4, 0.5, 0.6),
        (1, 1, 1)
    ]

    out_of_gamut_indices = []

    print("Analyzing ROMM RGB to sRGB conversion:\n")

    # --- 3. Process Each Color ---
    for i, rgb in enumerate(colors):
        color_num = i + 1
        romm_r, romm_g, romm_b = rgb

        # Delinearize the input ROMM RGB values
        linear_romm_vec = np.array([
            delinearize_romm_component(romm_r),
            delinearize_romm_component(romm_g),
            delinearize_romm_component(romm_b)
        ])

        # Transform to linear sRGB
        linear_srgb_vec = M_final_transform @ linear_romm_vec
        srgb_r, srgb_g, srgb_b = linear_srgb_vec

        # Check if the color is within the sRGB gamut [0, 1]
        is_in_gamut = (0 <= srgb_r <= 1) and (0 <= srgb_g <= 1) and (0 <= srgb_b <= 1)

        print(f"{color_num}) ROMM RGB({romm_r}, {romm_g}, {romm_b})")
        print(f"   -> Converts to linear sRGB({srgb_r:.4f}, {srgb_g:.4f}, {srgb_b:.4f})")
        
        if not is_in_gamut:
            out_of_gamut_indices.append(color_num)
            print("   -> Result: Out of sRGB gamut.\n")
        else:
            print("   -> Result: In sRGB gamut.\n")

    # --- 4. Final Answer ---
    print("-" * 30)
    if not out_of_gamut_indices:
        answer = "none cannot"
    else:
        answer = ", ".join(map(str, sorted(out_of_gamut_indices)))
    
    print(f"The numbers of the colors that cannot be represented are: {answer}")
    print(f"<<<{answer}>>>")

solve()