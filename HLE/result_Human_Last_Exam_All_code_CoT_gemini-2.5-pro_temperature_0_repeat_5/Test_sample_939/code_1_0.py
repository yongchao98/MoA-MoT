import numpy as np

def solve_color_gamut():
    """
    Determines which of the given ROMM RGB color values fall outside the sRGB gamut.
    """
    # Define the list of ROMM RGB colors to test (as a dictionary for easy labeling).
    colors = {
        1: (0, 0, 1),
        2: (0, 1, 0),
        3: (0, 0.5, 0.6),
        4: (0.4, 0.5, 0.6),
        5: (1, 1, 1)
    }

    # --- Transformation Matrices ---

    # Matrix to convert linear ROMM RGB (D50) to CIE XYZ (D50)
    M_ROMM_to_XYZ = np.array([
        [0.7976749, 0.1351917, 0.0313534],
        [0.2880402, 0.7118741, 0.0000857],
        [0.0000000, 0.0000000, 0.8252100]
    ])

    # Bradford Chromatic Adaptation Matrix from D50 to D65
    M_CAT_D50_to_D65 = np.array([
        [ 0.9555766, -0.0230393, 0.0631636],
        [-0.0282895,  1.0099416, 0.0210077],
        [ 0.0122982, -0.0204830, 1.3299098]
    ])

    # Matrix to convert CIE XYZ (D65) to linear sRGB (D65)
    M_XYZ_to_sRGB = np.array([
        [ 3.2404542, -1.5371385, -0.4985314],
        [-0.9692660,  1.8760108,  0.0415560],
        [ 0.0556434, -0.2040259,  1.0572252]
    ])

    out_of_gamut_indices = []

    # Iterate through each color and perform the conversion and check
    for index, romm_rgb in sorted(colors.items()):
        # Step 1: Apply ROMM RGB inverse gamma correction (gamma 1.8)
        romm_linear = np.array([c ** 1.8 if c > 0 else 0 for c in romm_rgb])

        # Step 2: Convert from linear ROMM (D50) to XYZ (D50)
        xyz_d50 = np.dot(M_ROMM_to_XYZ, romm_linear)

        # Step 3: Chromatically adapt from D50 to D65 white point
        xyz_d65 = np.dot(M_CAT_D50_to_D65, xyz_d50)

        # Step 4: Convert from XYZ (D65) to linear sRGB (D65)
        srgb_linear = np.dot(M_XYZ_to_sRGB, xyz_d65)
        
        r_srgb, g_srgb, b_srgb = srgb_linear
        print(f"Checking color {index}) ROMM RGB{romm_rgb}:")
        print(f"The final numbers in the conversion equation result in linear sRGB({r_srgb:.4f}, {g_srgb:.4f}, {b_srgb:.4f})")

        # Step 5: Check if the color is within the sRGB gamut [0, 1]
        # A small tolerance is used to account for floating-point inaccuracies.
        tolerance = 1e-9
        if np.any(srgb_linear < -tolerance) or np.any(srgb_linear > 1 + tolerance):
            out_of_gamut_indices.append(index)
            print("Result: Cannot be represented (out of sRGB gamut).\n")
        else:
            print("Result: Can be represented (in sRGB gamut).\n")

    # Format the final answer as requested
    if not out_of_gamut_indices:
        final_answer = "none cannot"
    else:
        final_answer = ", ".join(map(str, sorted(out_of_gamut_indices)))

    print(f"The numbers of the colors that cannot be represented are: {final_answer}")
    print(f"<<<{final_answer}>>>")

# Execute the function
solve_color_gamut()