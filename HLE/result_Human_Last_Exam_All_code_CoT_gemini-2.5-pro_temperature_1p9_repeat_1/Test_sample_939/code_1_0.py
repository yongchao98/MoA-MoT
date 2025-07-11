import numpy as np

def solve():
    """
    Checks which ROMM RGB colors are outside the sRGB gamut.
    """

    # Define the colors to be tested
    colors = {
        1: np.array([0, 0, 1]),
        2: np.array([0, 1, 0]),
        3: np.array([0, 0.5, 0.6]),
        4: np.array([0.4, 0.5, 0.6]),
        5: np.array([1, 1, 1]),
    }

    # 1. ProPhoto RGB (ROMM) primaries with D50 white point to XYZ D50
    M_romm_to_xyz_d50 = np.array([
        [0.7976749, 0.1351917, 0.0313534],
        [0.2880402, 0.7118741, 0.0000857],
        [0.0000000, 0.0000000, 0.8252100]
    ])

    # 2. Bradford Chromatic Adaptation Matrix from D50 to D65
    M_cat_d50_to_d65 = np.array([
        [ 0.9555766, -0.0230393, 0.0631636],
        [-0.0282895,  1.0099416, 0.0210077],
        [ 0.0122982, -0.0204830, 1.3299098]
    ])
    
    # 3. XYZ D65 to linear sRGB D65
    M_xyz_to_srgb_d65 = np.array([
        [ 3.24096994, -1.53738318, -0.49861076],
        [-0.96924364,  1.8759675,   0.04155506],
        [ 0.05563001, -0.20397696,  1.05697151]
    ])

    # Combine matrices for a single transformation from linear ROMM to linear sRGB
    M_final = M_xyz_to_srgb_d65 @ M_cat_d50_to_d65 @ M_romm_to_xyz_d50

    romm_gamma = 1.8
    out_of_gamut_indices = []

    print("Checking which ROMM RGB values can be represented in sRGB...")
    for index, romm_rgb in sorted(colors.items()):
        
        # Step 1: Linearize ROMM RGB (Gamma decoding)
        linear_romm = np.power(romm_rgb, romm_gamma)
        
        # Step 2: Convert to linear sRGB using the combined matrix
        linear_srgb = M_final @ linear_romm
        
        r, g, b = linear_srgb
        
        # Step 3: Check if the color is in the sRGB gamut [0, 1]
        is_in_gamut = (r >= 0 and r <= 1 and g >= 0 and g <= 1 and b >= 0 and b <= 1)
        
        result_str = "Can be represented" if is_in_gamut else "Cannot be represented"
        
        # Output the equation and result
        # Note: The 'equation' shows the final linear sRGB values which determine representability.
        print(f"{index}) ROMM RGB{tuple(romm_rgb)} -> linear sRGB({r:.4f}, {g:.4f}, {b:.4f}) -> {result_str}")

        if not is_in_gamut:
            out_of_gamut_indices.append(index)

    # Phrase the final answer as requested
    if not out_of_gamut_indices:
        final_answer = "none cannot"
    else:
        final_answer = ", ".join(map(str, sorted(out_of_gamut_indices)))
        
    print("\nThe numbers of the colors that cannot be represented are:")
    print(final_answer)
    print(f"<<<{final_answer}>>>")


solve()