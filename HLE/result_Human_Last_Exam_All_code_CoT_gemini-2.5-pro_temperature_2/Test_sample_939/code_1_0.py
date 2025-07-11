import numpy as np

def solve_color_gamut():
    """
    Analyzes which ROMM RGB colors cannot be represented in sRGB.
    """

    # Conversion matrix from linear ROMM RGB (D50) to CIE XYZ (D50)
    ROMM_TO_XYZ_D50 = np.array([
        [0.7976749, 0.1351917, 0.0313534],
        [0.2880402, 0.7118741, 0.0000857],
        [0.0000000, 0.0000000, 0.8252100]
    ])

    # Conversion matrix from CIE XYZ (D50) to linear sRGB (D65)
    # This matrix includes the Bradford chromatic adaptation transform.
    XYZ_D50_TO_SRGB_D65 = np.array([
        [ 3.1338561, -1.6168667, -0.4906146],
        [-0.9787684,  1.9161415,  0.0334540],
        [ 0.0719453, -0.2289914,  1.4052427]
    ])
    
    # Input ROMM RGB colors
    colors = {
        1: (0, 0, 1),
        2: (0, 1, 0),
        3: (0, 0.5, 0.6),
        4: (0.4, 0.5, 0.6),
        5: (1, 1, 1)
    }

    cannot_represent_indices = []

    print("Checking which ROMM RGB colors can be represented in sRGB...")

    for i, rgb_romm in colors.items():
        # --- The Conversion Equation Starts Here ---
        
        # Step 1: Linearize ROMM RGB (Inverse Companding)
        # ROMM uses a gamma of 1.8. Values are raised to the power of 1.8.
        # Note: A more precise definition has a linear segment near black, but
        # for these specific input values, a simple power function is sufficient.
        romm_r, romm_g, romm_b = rgb_romm
        
        # Applying the formula: V_lin = V_romm ^ 1.8
        linear_romm = np.power(np.array(rgb_romm), 1.8)
        
        # Step 2: Convert to XYZ
        # Applying the formula: [X,Y,Z] = Matrix * [R,G,B]_lin_romm
        xyz = ROMM_TO_XYZ_D50 @ linear_romm

        # Step 3: Convert to linear sRGB
        # Applying the formula: [R,G,B]_lin_srgb = Matrix * [X,Y,Z]
        srgb_linear = XYZ_D50_TO_SRGB_D65 @ xyz

        # --- The Conversion Equation Ends Here ---
        
        srgb_r, srgb_g, srgb_b = srgb_linear

        print(f"\n{i}) For ROMM RGB({romm_r}, {romm_g}, {romm_b}):")
        print(f"   The resulting linear sRGB values are: "
              f"R={srgb_r:.3f}, G={srgb_g:.3f}, B={srgb_b:.3f}")

        # Step 4: Check if the color is in sRGB gamut
        # The color is out of gamut if any component is < 0 or > 1.
        if np.any(srgb_linear < 0) or np.any(srgb_linear > 1):
            cannot_represent_indices.append(i)
            print("   Conclusion: This color IS NOT in the sRGB gamut.")
        else:
            print("   Conclusion: This color IS in the sRGB gamut.")

    print("\n---" * 10)
    if not cannot_represent_indices:
        result = "none cannot"
    else:
        # Format the final list as requested
        result = ", ".join(map(str, sorted(cannot_represent_indices)))

    print(f"The numbers of the colors that cannot be represented by an sRGB hex code are:")
    print(result)
    print(f"\n<<< {result} >>>")


solve_color_gamut()