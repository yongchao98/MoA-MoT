import numpy as np

def solve():
    """
    Determines which ROMM RGB colors cannot be represented in the sRGB color space.
    """
    # Conversion matrix from Linear ROMM RGB (D50) to Linear sRGB (D65).
    # This single matrix is a combination of three steps:
    # 1. ROMM RGB to CIE XYZ D50
    # 2. Chromatic Adaptation Transform (CAT02) from D50 to D65
    # 3. CIE XYZ D65 to sRGB
    M_ROMM_to_SRGB = np.array([
        [ 1.966993, -0.635914, -0.265103],
        [-0.231435,  1.163359,  0.076413],
        [-0.010867,  0.031780,  0.718878]
    ])

    colors = {
        1: np.array([0, 0, 1]),
        2: np.array([0, 1, 0]),
        3: np.array([0, 0.5, 0.6]),
        4: np.array([0.4, 0.5, 0.6]),
        5: np.array([1, 1, 1])
    }

    out_of_gamut_indices = []

    print("Checking which ROMM RGB colors are outside the sRGB gamut...")
    print("-" * 60)

    for i in sorted(colors.keys()):
        romm_rgb = colors[i]
        print(f"Checking color {i}: ROMM RGB({romm_rgb[0]}, {romm_rgb[1]}, {romm_rgb[2]})")

        # Step 1: Linearize the non-linear ROMM RGB values (inverse gamma correction)
        # For ROMM RGB, the gamma is 1.8.
        linear_romm_rgb = np.power(romm_rgb, 1.8)
        print(f"  - Linearized ROMM values: [{linear_romm_rgb[0]:.4f}, {linear_romm_rgb[1]:.4f}, {linear_romm_rgb[2]:.4f}]")

        # Step 2: Convert from linear ROMM RGB to linear sRGB using the matrix
        linear_srgb = M_ROMM_to_SRGB @ linear_romm_rgb
        print(f"  - Converted to Linear sRGB: [{linear_srgb[0]:.4f}, {linear_srgb[1]:.4f}, {linear_srgb[2]:.4f}]")

        # Step 3: Check if the linear sRGB values are in the [0, 1] gamut
        # We use a small tolerance for floating point inaccuracies.
        is_in_gamut = np.all(linear_srgb >= -1e-6) and np.all(linear_srgb <= 1 + 1e-6)
        
        reasons = []
        if linear_srgb[0] < 0: reasons.append(f"R={linear_srgb[0]:.2f} < 0")
        if linear_srgb[0] > 1: reasons.append(f"R={linear_srgb[0]:.2f} > 1")
        if linear_srgb[1] < 0: reasons.append(f"G={linear_srgb[1]:.2f} < 0")
        if linear_srgb[1] > 1: reasons.append(f"G={linear_srgb[1]:.2f} > 1")
        if linear_srgb[2] < 0: reasons.append(f"B={linear_srgb[2]:.2f} < 0")
        if linear_srgb[2] > 1: reasons.append(f"B={linear_srgb[2]:.2f} > 1")

        if is_in_gamut:
            print("  - Result: In sRGB gamut.")
        else:
            print(f"  - Result: Out of sRGB gamut. Reason(s): {', '.join(reasons)}")
            out_of_gamut_indices.append(str(i))
        
        print("-" * 60)

    if not out_of_gamut_indices:
        final_answer = "none cannot"
    else:
        final_answer = ", ".join(out_of_gamut_indices)

    print(f"The numbers of the colors that cannot be represented are: {final_answer}")
    print(f"<<<{final_answer}>>>")

solve()