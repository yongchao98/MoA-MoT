import numpy as np

def convert_and_check_gamut():
    """
    This function converts a list of ROMM RGB colors to sRGB and checks
    if they are within the sRGB gamut.
    """
    # Define transformation matrices for color space conversion
    # ROMM RGB (D50) to XYZ (D50)
    M_ROMM_TO_XYZ = np.array([
        [0.7976749, 0.1351917, 0.0313534],
        [0.2880402, 0.7118741, 0.0000857],
        [0.0000000, 0.0000000, 0.8252100]
    ])
    # Chromatic Adaptation: D50 to D65 (Bradford method)
    M_D50_TO_D65 = np.array([
        [0.9555766, -0.0230393, 0.0631636],
        [-0.0282895,  1.0099416, 0.0210077],
        [0.0122982, -0.0204830, 1.3299098]
    ])
    # XYZ (D65) to linear sRGB
    M_XYZ_TO_SRGB_LINEAR = np.array([
        [ 3.2404542, -1.5371385, -0.4985314],
        [-0.9692660,  1.8760108,  0.0415560],
        [ 0.0556434, -0.2040259,  1.0572252]
    ])

    # ROMM RGB uses a gamma of 1.8. The function to linearize it:
    def romm_to_linear(c_romm):
        # A small linear segment exists near black, but for the given values,
        # a simple power function is sufficient.
        if c_romm < (16 * 0.001953):
            return c_romm / 16.0
        else:
            return c_romm ** 1.8

    # The list of ROMM RGB colors to test
    colors = {
        1: (0, 0, 1),
        2: (0, 1, 0),
        3: (0, 0.5, 0.6),
        4: (0.4, 0.5, 0.6),
        5: (1, 1, 1)
    }

    out_of_gamut_indices = []

    print("Checking which ROMM RGB colors are outside the sRGB gamut:")
    for i, rgb_romm in colors.items():
        # Step 1: Linearize ROMM RGB
        rgb_romm_linear = np.array([romm_to_linear(c) for c in rgb_romm])
        
        # Step 2 & 3 & 4: Convert to linear sRGB
        xyz_d50 = np.dot(M_ROMM_TO_XYZ, rgb_romm_linear)
        xyz_d65 = np.dot(M_D50_TO_D65, xyz_d50)
        srgb_linear = np.dot(M_XYZ_TO_SRGB_LINEAR, xyz_d65)

        # Step 5: Check gamut. A color is out of gamut if any component is < 0 or > 1.
        # A small tolerance is used for floating-point inaccuracies.
        tolerance = 1e-4
        if np.any(srgb_linear < -tolerance) or np.any(srgb_linear > 1.0 + tolerance):
            is_in_gamut = False
            out_of_gamut_indices.append(i)
        else:
            is_in_gamut = True
        
        status = "IN sRGB gamut" if is_in_gamut else "OUT of sRGB gamut"
        print(f"{i}) ROMM RGB {rgb_romm} -> Linear sRGB ({srgb_linear[0]:.3f}, {srgb_linear[1]:.3f}, {srgb_linear[2]:.3f}) -> {status}")

    if not out_of_gamut_indices:
        result = "none cannot"
    else:
        out_of_gamut_indices.sort()
        result = ", ".join(map(str, out_of_gamut_indices))

    print("\n" + "="*30)
    print("Final Answer")
    print("="*30)
    print(f"{result}")

convert_and_check_gamut()