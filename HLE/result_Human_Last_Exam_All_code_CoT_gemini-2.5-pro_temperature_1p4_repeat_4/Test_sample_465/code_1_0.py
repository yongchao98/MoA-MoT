import numpy as np
import numpy.polynomial.polynomial as poly

def demonstrate_chromatic_roots():
    """
    This function demonstrates properties of chromatic roots using numpy
    to find the roots of chromatic polynomials for specific graphs.
    """
    # --- Statement B: Chromatic roots may not be real (Example: C4 graph) ---
    print("--- Demonstrating Statement B (complex roots) with C4 Graph ---")
    # The chromatic polynomial of C4 is P(k) = k^4 - 4k^3 + 6k^2 - 3k.
    # The equation to solve is P(k) = 0.
    c4_coeffs = [1, -4, 6, -3, 0]
    print("The equation for C4 chromatic roots is:")
    print("1*k^4 - 4*k^3 + 6*k^2 - 3*k + 0 = 0")
    
    c4_roots = np.roots(c4_coeffs)
    print("\nChromatic roots of C4:")
    for root in c4_roots:
        # Using a small tolerance to check for purely real/imaginary parts
        if np.isclose(root.imag, 0):
            print(f"  {root.real:.4f} (Real)")
        else:
            print(f"  {root.real:.4f} + {root.imag:.4f}j (Complex)")
    print("As shown, C4 has complex chromatic roots, confirming statement B.\n")

    # --- Statement D: Real chromatic roots may be non-integers (Example: Petersen Graph) ---
    print("--- Demonstrating Statement D (non-integer real roots) with Petersen Graph ---")
    # The chromatic polynomial of the Petersen graph is P(k) = k(k-1)(k-2) * Q(k),
    # where Q(k) = k^7 - 12k^6 + 67k^5 - 230k^4 + 529k^3 - 814k^2 + 775k - 352.
    # We find the roots of Q(k) to find the non-trivial roots.
    petersen_factor_coeffs = [1, -12, 67, -230, 529, -814, 775, -352]
    print("The non-trivial roots of the Petersen graph's chromatic polynomial solve Q(k) = 0:")
    print("1*k^7 - 12*k^6 + 67*k^5 - 230*k^4 + 529*k^3 - 814*k^2 + 775*k - 352 = 0")
    
    petersen_roots = np.roots(petersen_factor_coeffs)
    print("\nNon-trivial chromatic roots of the Petersen Graph:")
    for root in petersen_roots:
        if np.isclose(root.imag, 0):
            # Check if the root is an integer
            if np.isclose(root.real, round(root.real)):
                 print(f"  {root.real:.4f} (Real Integer)")
            else:
                 print(f"  {root.real:.4f} (Real Non-Integer)")
        else:
            print(f"  {root.real:.4f} + {root.imag:.4f}j (Complex)")
    print("As shown, the Petersen graph has a real, non-integer root (~2.8536), confirming statement D.")

demonstrate_chromatic_roots()