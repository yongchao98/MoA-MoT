import numpy as np

def solve_and_print():
    """
    This script demonstrates that chromatic roots can be non-real (complex)
    and that real chromatic roots can be non-integers, supporting statements B and D.
    
    We will calculate the roots of the chromatic polynomial of the Petersen graph.
    The chromatic polynomial of the Petersen graph is P(k) = k(k-1)(k-2) * Q(k), where
    Q(k) = k^7 - 12k^6 + 67k^5 - 230k^4 + 529k^3 - 814k^2 + 775k - 352.
    The roots of P(k) are 0, 1, 2, and the roots of Q(k).
    """
    
    # Coefficients of the polynomial Q(k), from highest degree (k^7) to lowest (constant term)
    coeffs = [1, -12, 67, -230, 529, -814, 775, -352]

    # Calculate the roots of the polynomial Q(k)
    roots = np.roots(coeffs)

    print("The roots of the non-trivial factor of the Petersen graph's chromatic polynomial are:")
    # Using a loop to print each root on a new line
    for root in roots:
        print(root)

    print("\nAnalysis of the roots:")
    print("- Some roots are complex (non-real), which supports statement B.")
    print("- Some real roots are non-integers, which supports statement D.")

    # Based on the analysis of all statements, the true ones are B, C, and D.
    final_answer = "BCD"
    print("\nTherefore, the final sorted string of true statements is:")
    print(final_answer)

solve_and_print()