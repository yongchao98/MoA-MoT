import numpy as np

def solve():
    """
    Determines which ROMM RGB colors are outside the sRGB gamut.
    """
    # This matrix converts from linear ROMM RGB to linear sRGB,
    # including chromatic adaptation from D50 to D65.
    # It is derived from the standard transformations:
    # M_sRGB_from_XYZ @ M_CAT_D50_D65 @ M_XYZ_from_ROMM
    M_romm_to_srgb = np.array([
        [ 1.34579896, -0.2555801,  -0.05110593],
        [-0.54462249,  1.50823613, -0.0040523 ],
        [-0.00000571, -0.00406858,  0.69183244]
    ])

    # List of the ROMM RGB colors to check
    colors = [
        np.array([0, 0, 1]),      # 1
        np.array([0, 1, 0]),      # 2
        np.array([0, 0.5, 0.6]),  # 3
        np.array([0.4, 0.5, 0.6]),# 4
        np.array([1, 1, 1])       # 5
    ]

    out_of_gamut_indices = []

    print("Checking which ROMM RGB colors are outside the sRGB gamut...")
    print("-" * 60)

    for i, romm_nl in enumerate(colors):
        # 1. Linearize the non-linear ROMM RGB value (gamma 1.8)
        romm_lin = np.power(romm_nl, 1.8)

        # 2. Convert from linear ROMM RGB to linear sRGB
        srgb_lin = M_romm_to_srgb @ romm_lin

        # 3. Check if the color is within the sRGB gamut [0, 1]
        # A small tolerance is used for floating point comparisons
        tolerance = 1e-9
        is_in_gamut = np.all(srgb_lin >= -tolerance) and np.all(srgb_lin <= 1 + tolerance)

        print(f"Color {i+1}: ROMM RGB({romm_nl[0]}, {romm_nl[1]}, {romm_nl[2]})")
        print(f"  -> Linear sRGB({srgb_lin[0]:.4f}, {srgb_lin[1]:.4f}, {srgb_lin[2]:.4f})")
        
        if not is_in_gamut:
            out_of_gamut_indices.append(str(i + 1))
            print("  -> Result: Cannot be represented (Out of sRGB gamut)")
        else:
            print("  -> Result: Can be represented (In sRGB gamut)")
        print("-" * 60)

    if not out_of_gamut_indices:
        result = "none cannot"
    else:
        result = ", ".join(out_of_gamut_indices)

    print(f"\nThe numbers of the colors that cannot be represented are:\n{result}")
    
    # Final answer in the required format
    print(f"\n<<<{result}>>>")

solve()