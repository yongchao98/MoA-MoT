import numpy as np
import itertools

def analyze_stbc_diversity():
    """
    Analyzes the diversity order of three STBCs by checking the rank criterion.
    """
    # Use 4-QAM constellation: {1+j, 1-j, -1+j, -1-j}
    qam_symbols = [
        1 + 1j, 1 - 1j, -1 + 1j, -1 - 1j
    ]
    
    # Generate all possible difference symbols delta_x = x_i - x_j
    # We are interested in non-zero differences
    diff_symbols = set()
    for s1 in qam_symbols:
        for s2 in qam_symbols:
            diff_symbols.add(s1 - s2)

    # All possible pairs of (delta_x1, delta_x2) where at least one is non-zero
    diff_pairs = list(itertools.product(diff_symbols, repeat=2))
    # Filter out the (0,0) pair
    diff_pairs = [p for p in diff_pairs if p[0] != 0 or p[1] != 0]

    # --- Analysis for Code Sa ---
    print("--- Analysis for Code Sa: S = [[x1, x2], [x2, x1]] ---")
    found_singular_a = False
    for dx1, dx2 in diff_pairs:
        # det(Delta_S_a) = (dx1)^2 - (dx2)^2
        det_a = dx1**2 - dx2**2
        if np.isclose(det_a, 0):
            print("Found a singular difference matrix for Sa:")
            print(f"  delta_x1 = {dx1}")
            print(f"  delta_x2 = {dx2}")
            print(f"  Determinant = ({dx1})^2 - ({dx2})^2 = {dx1**2} - {dx2**2} = {det_a}")
            found_singular_a = True
            break
    if found_singular_a:
        print("Diversity order for Sa is 1.\n")
    else:
        print("No singular difference matrix found. Diversity order for Sa is 2.\n")

    # --- Analysis for Code Sb ---
    print("--- Analysis for Code Sb: S = [[x1, x2], [x2, x1*]] ---")
    found_singular_b = False
    for dx1, dx2 in diff_pairs:
        # det(Delta_S_b) = |dx1|^2 - (dx2)^2
        det_b = abs(dx1)**2 - dx2**2
        if np.isclose(det_b, 0):
            print("Found a singular difference matrix for Sb:")
            print(f"  delta_x1 = {dx1}")
            print(f"  delta_x2 = {dx2}")
            print(f"  Determinant = |{dx1}|^2 - ({dx2})^2 = {abs(dx1)**2:.1f} - {dx2**2} = {det_b}")
            found_singular_b = True
            break
    if found_singular_b:
        print("Diversity order for Sb is 1.\n")
    else:
        print("No singular difference matrix found. Diversity order for Sb is 2.\n")
        
    # --- Analysis for Code Sc ---
    print("--- Analysis for Code Sc: S = [[-x1*, x2], [-x2*, -x1]] ---")
    found_singular_c = False
    for dx1, dx2 in diff_pairs:
        # det(Delta_S_c) = |dx1|^2 + |dx2|^2
        det_c = abs(dx1)**2 + abs(dx2)**2
        if np.isclose(det_c, 0):
            print("Found a singular difference matrix for Sc:")
            print(f"  delta_x1 = {dx1}")
            print(f"  delta_x2 = {dx2}")
            print(f"  Determinant = |{dx1}|^2 + |{dx2}|^2 = {abs(dx1)**2:.1f} + {abs(dx2)**2:.1f} = {det_c}")
            found_singular_c = True
            break
            
    if found_singular_c:
        print("This should not happen for non-zero delta vectors.")
        print("Diversity order for Sc is 1.\n")
    else:
        # We can also show a sample calculation for a non-zero case
        dx1, dx2 = 2+0j, 2j
        det_c_sample = abs(dx1)**2 + abs(dx2)**2
        print(f"No singular difference matrix found for Sc over {len(diff_pairs)} pairs.")
        print(f"The determinant |delta_x1|^2 + |delta_x2|^2 is always positive.")
        print("Example calculation:")
        print(f"  delta_x1 = {dx1}")
        print(f"  delta_x2 = {dx2}")
        print(f"  Determinant = |{dx1}|^2 + |{dx2}|^2 = {abs(dx1)**2:.1f} + {abs(dx2)**2:.1f} = {det_c_sample:.1f}")
        print("Diversity order for Sc is 2.\n")
        
    # --- Final Answers ---
    print("="*30)
    print("Final Answers:")
    print("(a) What is the diversity order for each code?")
    print("    - Diversity order of Sa is 1.")
    print("    - Diversity order of Sb is 1.")
    print("    - Diversity order of Sc is 2.")
    print("\n(b) Which code provides the maximum directivity (diversity)?")
    print("    Code Sc provides the maximum diversity order of 2.")
    print("="*30)

if __name__ == '__main__':
    analyze_stbc_diversity()
