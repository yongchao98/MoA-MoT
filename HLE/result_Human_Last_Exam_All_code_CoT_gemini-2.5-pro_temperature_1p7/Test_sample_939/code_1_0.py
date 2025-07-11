import numpy as np

def romm_to_srgb_check():
    """
    Checks which ROMM RGB colors cannot be represented in sRGB.

    This function converts a list of ROMM RGB colors to the sRGB color space
    and checks if they fall within the sRGB gamut [0, 1].
    """
    
    # Standard transformation matrix from linear ROMM RGB (D50) to linear sRGB (D65)
    # using the Bradford chromatic adaptation method.
    M_ROMM_TO_SRGB = np.array([
        [1.2536831, -0.0768913, -0.1537234],
        [-0.2435749,  1.1784381,  0.0652576],
        [-0.0093836, -0.2941544,  1.3094896]
    ])

    colors = {
        1: (0, 0, 1),
        2: (0, 1, 0),
        3: (0, 0.5, 0.6),
        4: (0.4, 0.5, 0.6),
        5: (1, 1, 1)
    }
    
    cannot_represent = []

    print("Checking which ROMM RGB values are outside the sRGB gamut...")
    print("-" * 60)

    for number, rgb_nl in colors.items():
        # Step 1: Convert non-linear ROMM RGB to linear values
        # The transfer function for ROMM RGB (inverse gamma) is applied.
        romm_lin = []
        for val_nl in rgb_nl:
            if val_nl >= 0.03125:
                romm_lin.append(val_nl ** 1.8)
            else:
                romm_lin.append(val_nl / 16.0)
        
        romm_lin_vec = np.array(romm_lin)

        # Step 2: Apply the transformation matrix
        srgb_lin_vec = M_ROMM_TO_SRGB @ romm_lin_vec

        # Step 3: Check if the resulting sRGB color is in gamut [0, 1]
        # We use a small tolerance for the upper bound to account for
        # floating point inaccuracies, especially for the white point (1,1,1).
        # Any negative value is definitively out of gamut.
        is_in_gamut = all(component >= -1e-5 and component <= 1.0 + 1e-2 for component in srgb_lin_vec)
        
        # Format for clear output
        romm_str = f"ROMM RGB({rgb_nl[0]}, {rgb_nl[1]}, {rgb_nl[2]})"
        srgb_str = f"sRGB(R={srgb_lin_vec[0]:.3f}, G={srgb_lin_vec[1]:.3f}, B={srgb_lin_vec[2]:.3f})"
        
        if is_in_gamut:
            print(f"{number}) {romm_str} -> {srgb_str} -> Can be represented.")
        else:
            print(f"{number}) {romm_str} -> {srgb_str} -> CANNOT be represented.")
            cannot_represent.append(number)
    
    print("-" * 60)
    
    if not cannot_represent:
        result = "none cannot"
    else:
        result = ", ".join(map(str, sorted(cannot_represent)))
        
    print("Numbers of the colors that cannot be represented:")
    print(result)
    
    # Final answer in the required format
    print(f"<<<{result}>>>")

if __name__ == "__main__":
    romm_to_srgb_check()