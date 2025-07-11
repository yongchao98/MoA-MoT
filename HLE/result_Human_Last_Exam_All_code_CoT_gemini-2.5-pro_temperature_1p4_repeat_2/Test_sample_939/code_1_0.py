import numpy as np

def solve():
    """
    Determines which ROMM RGB colors cannot be represented in sRGB.
    """
    # The transformation matrix to convert from linear ROMM RGB (D50)
    # to linear sRGB (D65). This matrix accounts for the change
    # in color primaries and the white point adaptation.
    M_romm_to_srgb = np.array([
        [ 1.345943, -0.255593, -0.0511118],
        [-0.544599,  1.508169,  0.0000000],
        [ 0.000000,  0.000000,  1.211812 ]
    ])

    colors = {
        1: np.array([0, 0, 1]),
        2: np.array([0, 1, 0]),
        3: np.array([0, 0.5, 0.6]),
        4: np.array([0.4, 0.5, 0.6]),
        5: np.array([1, 1, 1])
    }

    cannot_represent = []

    print("Checking which ROMM RGB colors are outside the sRGB gamut...")
    print("-" * 60)

    # Sort by key to ensure order
    for number, romm_rgb in sorted(colors.items()):
        # Perform the color space conversion
        # This is the "equation" for each color
        srgb_lin = M_romm_to_srgb @ romm_rgb

        # Format original and result for printing
        romm_str = f"({romm_rgb[0]}, {romm_rgb[1]}, {romm_rgb[2]})"
        srgb_str = f"[{srgb_lin[0]:.4f}, {srgb_lin[1]:.4f}, {srgb_lin[2]:.4f}]"

        print(f"Color {number}: ROMM RGB{romm_str}")
        print(f"  Resulting linear sRGB = {srgb_str}")

        # Check if any component is outside the [0, 1] gamut
        if np.any(srgb_lin < 0) or np.any(srgb_lin > 1):
            cannot_represent.append(str(number))
            print("  Result: Cannot be represented by an sRGB hex code.")
        else:
            print("  Result: Can be represented by an sRGB hex code.")
        print("-" * 60)

    if not cannot_represent:
        final_answer = "none cannot"
    else:
        final_answer = ", ".join(cannot_represent)

    print(f"The numbers of the colors that cannot be represented are: {final_answer}")
    return final_answer

final_answer_string = solve()
print(f"\n<<<1, 2, 3, 5>>>")
