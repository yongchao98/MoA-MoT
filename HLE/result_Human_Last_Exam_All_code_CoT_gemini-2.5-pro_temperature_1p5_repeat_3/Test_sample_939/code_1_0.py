import numpy as np

def solve():
    """
    Checks which ROMM RGB colors are outside the sRGB gamut.
    """

    # This matrix converts linear ROMM RGB (D50 white point)
    # to linear sRGB (D65 white point).
    # It is derived from the standard ProPhoto-to-XYZ, Bradford Chromatic
    # Adaptation (D50 to D65), and XYZ-to-sRGB matrices.
    ROMM_TO_SRGB_MATRIX = np.array([
        [ 1.84180497, -0.53320078, -0.26449195],
        [-0.98501258,  1.91336441,  0.03328211],
        [ 0.05755191, -0.12999547,  1.1399638 ]
    ])

    colors = {
        1: np.array([0, 0, 1]),
        2: np.array([0, 1, 0]),
        3: np.array([0, 0.5, 0.6]),
        4: np.array([0.4, 0.5, 0.6]),
        5: np.array([1, 1, 1])
    }

    cannot_represent = []

    print("Analyzing ROMM RGB to sRGB conversion:\n")

    for num, romm_rgb in sorted(colors.items()):
        # Step 1: Linearize the ROMM RGB values by applying inverse gamma (1.8)
        # Note: ROMM RGB standard has a linear segment for dark colors, but for the
        # values in this problem, a simple power function is sufficient.
        linear_romm_rgb = np.power(romm_rgb, 1.8)

        # Step 2: Convert from linear ROMM RGB to linear sRGB using the matrix
        linear_srgb = ROMM_TO_SRGB_MATRIX.dot(linear_romm_rgb)

        # Step 3: Check if the resulting sRGB values are in gamut [0, 1]
        is_in_gamut = np.all(linear_srgb >= 0) and np.all(linear_srgb <= 1)
        
        # --- Output the calculation details ---
        print(f"Color {num}: ROMM({romm_rgb[0]}, {romm_rgb[1]}, {romm_rgb[2]})")
        print(f"  Linearized ROMM: ({linear_romm_rgb[0]:.4f}, {linear_romm_rgb[1]:.4f}, {linear_romm_rgb[2]:.4f})")
        
        # Unpack values for equation clarity
        r_romm, g_romm, b_romm = linear_romm_rgb
        m = ROMM_TO_SRGB_MATRIX
        r_srgb, g_srgb, b_srgb = linear_srgb
        
        print(f"  Converting to linear sRGB:")
        print(f"    sR = ({m[0,0]:.2f} * {r_romm:.2f}) + ({m[0,1]:.2f} * {g_romm:.2f}) + ({m[0,2]:.2f} * {b_romm:.2f}) = {r_srgb:.4f}")
        print(f"    sG = ({m[1,0]:.2f} * {r_romm:.2f}) + ({m[1,1]:.2f} * {g_romm:.2f}) + ({m[1,2]:.2f} * {b_romm:.2f}) = {g_srgb:.4f}")
        print(f"    sB = ({m[2,0]:.2f} * {r_romm:.2f}) + ({m[2,1]:.2f} * {g_romm:.2f}) + ({m[2,2]:.2f} * {b_romm:.2f}) = {b_srgb:.4f}")

        if is_in_gamut:
            print("  Result: In sRGB gamut. CAN be represented.\n")
        else:
            print("  Result: OUT of sRGB gamut. CANNOT be represented.\n")
            cannot_represent.append(str(num))

    final_answer = ", ".join(cannot_represent)
    print(f"The numbers of the colors that cannot be represented are: {final_answer}")
    print(f"<<<{final_answer}>>>")

solve()