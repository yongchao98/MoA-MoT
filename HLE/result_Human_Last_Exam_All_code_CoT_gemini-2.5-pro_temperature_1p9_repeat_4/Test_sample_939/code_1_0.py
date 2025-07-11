import numpy as np

def solve():
    """
    Checks which ROMM RGB colors from a given list are outside the sRGB gamut.
    """
    # List of ROMM RGB color values to check
    colors = [
        (0, 0, 1),      # 1
        (0, 1, 0),      # 2
        (0, 0.5, 0.6),  # 3
        (0.4, 0.5, 0.6),# 4
        (1, 1, 1)       # 5
    ]

    # --- Transformation Matrices ---

    # 1. ROMM RGB (ProPhoto) to XYZ (D50 illuminant)
    M_ROMM_to_XYZ = np.array([
        [0.7976749, 0.1351917, 0.0313534],
        [0.2880402, 0.7118741, 0.0000857],
        [0.0000000, 0.0000000, 0.8252100]
    ])

    # 2. Chromatic Adaptation Matrix from D50 to D65 (Bradford method)
    M_D50_to_D65 = np.array([
        [ 1.0478112,  0.0228866, -0.0501270],
        [ 0.0295424,  0.9904844, -0.0170491],
        [-0.0092345,  0.0150453,  0.7521316]
    ])

    # 3. XYZ (D65 illuminant) to linear sRGB
    M_XYZ_to_sRGB = np.array([
        [ 3.2404542, -1.5371385, -0.4985314],
        [-0.9692660,  1.8760108,  0.0415560],
        [ 0.0556434, -0.2040259,  1.0572252]
    ])

    cannot_represent = []

    print("Checking which ROMM RGB values are outside the sRGB gamut.")
    print("A color is out of gamut if its linear sRGB components are outside the [0, 1] range.\n")


    for i, color_romm in enumerate(colors):
        # Step 1: Linearize ROMM RGB. The given values are non-linear (gamma corrected).
        # For the given inputs, this is done by applying an exponent of 1.8.
        # This is a simplification; a full implementation has a linear segment for
        # very dark values, but none of our inputs (except 0) are in that segment.
        r_lin = color_romm[0] ** 1.8 if color_romm[0] > 0 else 0
        g_lin = color_romm[1] ** 1.8 if color_romm[1] > 0 else 0
        b_lin = color_romm[2] ** 1.8 if color_romm[2] > 0 else 0
        linear_romm = np.array([r_lin, g_lin, b_lin])

        # Step 2, 3, 4: Convert through XYZ to linear sRGB
        xyz_d50 = M_ROMM_to_XYZ @ linear_romm
        xyz_d65 = M_D50_to_D65 @ xyz_d50
        linear_srgb = M_XYZ_to_sRGB @ xyz_d65
        
        # Unpack for printing the final equation values
        R, G, B = linear_srgb[0], linear_srgb[1], linear_srgb[2]

        # Step 5: Check if the color is within the sRGB gamut
        # The condition for being in-gamut is 0 <= component <= 1 for all components.
        is_in_gamut = np.all((linear_srgb >= 0) & (linear_srgb <= 1))

        print(f"Color {i+1}: ROMM{color_romm}")
        print(f"  - Final linear sRGB values are R={R:.4f}, G={G:.4f}, B={B:.4f}")

        if not is_in_gamut:
            cannot_represent.append(i + 1)
            print(f"  - Status: CANNOT be represented in sRGB.")
        else:
            print(f"  - Status: Can be represented in sRGB.")
        print("-" * 30)


    # Final result
    if not cannot_represent:
        result = "none cannot"
    else:
        # Format the final list of numbers as requested
        result = ", ".join(map(str, sorted(cannot_represent)))

    print("\nNumbers of the colors that cannot be represented:")
    print(result)


solve()