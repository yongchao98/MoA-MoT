import numpy as np

def solve():
    """
    Determines which ROMM RGB colors cannot be represented in sRGB.
    """

    # List of ROMM RGB colors to test, with their corresponding number
    colors = {
        1: np.array([0, 0, 1]),
        2: np.array([0, 1, 0]),
        3: np.array([0, 0.5, 0.6]),
        4: np.array([0.4, 0.5, 0.6]),
        5: np.array([1, 1, 1]),
    }

    # This is the standard matrix for converting linear ROMM RGB (ProPhoto primaries, D50 white point)
    # to linear sRGB (sRGB primaries, D65 white point), including chromatic adaptation.
    M_ROMM_to_sRGB = np.array([
        [ 1.46200135, -0.18490947, -0.27436034],
        [-0.52181313,  1.44723812,  0.06786193],
        [ 0.03495449, -0.09689814,  1.2884838 ]
    ])

    cannot_represent = []

    print("Analyzing ROMM RGB to sRGB conversion for each color:")
    print("-" * 65)

    # Sort by the color number to ensure ascending order
    for i in sorted(colors.keys()):
        rgb_romm = colors[i]
        
        # Start of equation output for the current color
        print(f"Color {i}: RGB({rgb_romm[0]}, {rgb_romm[1]}, {rgb_romm[2]})")

        # Step 1: Linearize ROMM RGB values. Linear value = (ROMM value) ^ 1.8
        rgb_romm_linear = np.power(rgb_romm, 1.8)
        print(f"  1. Linearized ROMM Values:        ({rgb_romm_linear[0]:.4f}, {rgb_romm_linear[1]:.4f}, {rgb_romm_linear[2]:.4f})")

        # Step 2: Convert from linear ROMM RGB to linear sRGB using the matrix
        rgb_srgb_linear = M_ROMM_to_sRGB.dot(rgb_romm_linear)
        print(f"  2. Resulting Linear sRGB (R,G,B): ({rgb_srgb_linear[0]:.4f}, {rgb_srgb_linear[1]:.4f}, {rgb_srgb_linear[2]:.4f})")

        # Step 3: Check if the color is within the sRGB gamut [0, 1].
        # Add a small tolerance for floating point inaccuracies.
        is_in_gamut = np.all((rgb_srgb_linear >= -0.0001) & (rgb_srgb_linear <= 1.0001)) 
        
        reasons = []
        if rgb_srgb_linear[0] < 0.0: reasons.append("R < 0")
        if rgb_srgb_linear[0] > 1.0: reasons.append("R > 1")
        if rgb_srgb_linear[1] < 0.0: reasons.append("G < 0")
        if rgb_srgb_linear[1] > 1.0: reasons.append("G > 1")
        if rgb_srgb_linear[2] < 0.0: reasons.append("B < 0")
        if rgb_srgb_linear[2] > 1.0: reasons.append("B > 1")

        if is_in_gamut:
            print("  3. Result: Can be represented in sRGB.")
        else:
            print(f"  3. Result: Cannot be represented in sRGB ({', '.join(reasons)}).")
            cannot_represent.append(str(i))
        
        print("-" * 65)

    if not cannot_represent:
        final_answer = "none cannot"
    else:
        final_answer = ", ".join(cannot_represent)

    print(f"\nFinal Answer: The numbers of the colors that cannot be represented are {final_answer}")
    print(f"<<<{final_answer}>>>")

solve()