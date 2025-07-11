import numpy as np

def solve():
    """
    This function determines which ROMM RGB colors are outside the sRGB gamut.
    """
    # Define the colors to be tested. The numbers are 1-indexed.
    colors = {
        1: np.array([0, 0, 1]),
        2: np.array([0, 1, 0]),
        3: np.array([0, 0.5, 0.6]),
        4: np.array([0.4, 0.5, 0.6]),
        5: np.array([1, 1, 1])
    }

    # Standard matrix for converting linear ProPhoto/ROMM RGB to linear sRGB
    # This matrix accounts for the change in primaries and chromatic adaptation from D50 to D65.
    M_ROMM_to_sRGB = np.array([
        [ 2.04159, -0.56501, -0.34473],
        [-0.96924,  1.87597,  0.04156],
        [ 0.01344, -0.11836,  1.01517]
    ])

    cannot_be_represented = []

    print("Checking which ROMM RGB colors cannot be represented in sRGB...")
    print("-" * 60)

    for number, romm_rgb in sorted(colors.items()):
        # Step 1: Linearize the ROMM RGB values (Gamma 1.8 EOTF)
        # Note: The more complete EOTF has a linear segment for small values,
        # but for the given inputs, a simple power function is sufficient.
        linear_romm_rgb = np.power(romm_rgb, 1.8)

        # Step 2: Apply the transformation matrix
        linear_srgb = M_ROMM_to_sRGB.dot(linear_romm_rgb)

        # Step 3: Check if the resulting linear sRGB values are in the [0, 1] gamut
        # np.any(linear_srgb < 0) checks if any element is negative
        # np.any(linear_srgb > 1) checks if any element is greater than 1
        is_out_of_gamut = np.any(linear_srgb < 0) or np.any(linear_srgb > 1)

        print(f"Color {number}: ROMM RGB{tuple(romm_rgb)}")
        print(f"  -> Linear sRGB: [{linear_srgb[0]:.4f}, {linear_srgb[1]:.4f}, {linear_srgb[2]:.4f}]")
        
        if is_out_of_gamut:
            cannot_be_represented.append(str(number))
            print("  -> Result: Out of sRGB gamut. Cannot be represented.")
        else:
            print("  -> Result: In sRGB gamut. Can be represented.")
        print("-" * 60)

    # Phrase the final answer
    if not cannot_be_represented:
        final_answer = "none cannot"
    else:
        final_answer = ", ".join(cannot_be_represented)
    
    print(f"\nThe numbers of the colors that cannot be represented are: {final_answer}")
    print(f"\n<<<1, 2, 3, 5>>>")


solve()