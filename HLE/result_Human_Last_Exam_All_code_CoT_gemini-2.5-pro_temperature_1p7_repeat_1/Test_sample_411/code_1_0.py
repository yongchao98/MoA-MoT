import numpy as np

def analyze_stbc_diversity():
    """
    Analyzes and demonstrates the diversity order of three space-time codes.
    """
    print("--- Diversity Order Analysis ---")
    print("A diversity order of N=2 is achieved if the rank of the difference matrix dS is always 2.")
    print("If a case exists where rank(dS) = 1, the diversity order is 1.\n")

    # Define symbols from QPSK constellation
    s = 1 / np.sqrt(2)
    qpsk_symbols = [s * (1 + 1j), s * (-1 + 1j), s * (-1 - 1j), s * (1 - 1j)]

    # --- Code Sa ---
    print("Analysis of Code Sa: dS = [[dx1, dx2], [dx2, dx1]]")
    # We find an example where det(dSa) = (dx1)^2 - (dx2)^2 = 0, which occurs if dx1 = dx2.
    # Choose codewords (x1, x2) and (x1', x2') such that x1-x1' = x2-x2'.
    x1, x2 = qpsk_symbols[0], qpsk_symbols[0]
    x1_p, x2_p = qpsk_symbols[1], qpsk_symbols[1]
    dx1 = x1 - x1_p
    dx2 = x2 - x2_p
    dSa = np.array([[dx1, dx2], [dx2, dx1]])
    det_dSa = np.linalg.det(dSa)
    rank_dSa = np.linalg.matrix_rank(dSa, tol=1e-9)

    print(f"Let dx1 = x1-x1' = {dx1:.2f}")
    print(f"Let dx2 = x2-x2' = {dx2:.2f}")
    print(f"Equation: det(dSa) = (dx1)^2 - (dx2)^2 = ({dx1:.2f})^2 - ({dx2:.2f})^2 = {dx1**2 - dx2**2:.2f}")
    print(f"Matrix dSa:\n{dSa}")
    print(f"Result: det(dSa) is {det_dSa.real:.2f} and rank is {rank_dSa}. Diversity Order = 1.\n")

    # --- Code Sb ---
    print("Analysis of Code Sb: dS = [[dx1, dx2], [dx2, dx1*]]")
    # We find an example where det(dSb) = |dx1|^2 - (dx2)^2 = 0.
    dx1 = qpsk_symbols[0] - qpsk_symbols[3]  # results in j*sqrt(2)
    dx2 = qpsk_symbols[0] - qpsk_symbols[1]  # results in sqrt(2)
    dSb = np.array([[dx1, dx2], [dx2, np.conj(dx1)]])
    det_dSb = np.linalg.det(dSb)
    rank_dSb = np.linalg.matrix_rank(dSb, tol=1e-9)

    print(f"Let dx1 = {dx1:.2f} and dx2 = {dx2:.2f}")
    print(f"Equation: det(dSb) = |dx1|^2 - (dx2)^2 = {np.abs(dx1)**2:.1f} - ({dx2.real:.2f})^2 = {np.abs(dx1)**2 - dx2**2:.2f}")
    print(f"Matrix dSb:\n{dSb}")
    print(f"Result: det(dSb) is {det_dSb.real:.2f} and rank is {rank_dSb}. Diversity Order = 1.\n")
    
    # --- Code Sc ---
    print("Analysis of Code Sc: dS = [[-dx1*, dx2], [-dx2*, -dx1]]")
    # The determinant is det(dSc) = |dx1|^2 + |dx2|^2, which is never zero for dx1, dx2 != 0.
    # We use the same dx1 and dx2 from the Sb example to show the rank is 2.
    dSc = np.array([[-np.conj(dx1), dx2], [-np.conj(dx2), -dx1]])
    det_dSc = np.linalg.det(dSc)
    rank_dSc = np.linalg.matrix_rank(dSc, tol=1e-9)

    print(f"Using dx1 = {dx1:.2f} and dx2 = {dx2:.2f}")
    print(f"Equation: det(dSc) = |dx1|^2 + |dx2|^2 = {np.abs(dx1)**2:.1f} + {np.abs(dx2)**2:.1f} = {np.abs(dx1)**2 + np.abs(dx2)**2:.2f}")
    print(f"Matrix dSc:\n{dSc}")
    print(f"Result: det(dSc) is {det_dSc.real:.2f} and rank is {rank_dSc}. Diversity Order = 2.\n")
    
    print("--- Conclusion ---")
    print("(a) The diversity orders for Sa, Sb, and Sc are 1, 1, and 2, respectively.")
    print("(b) Code Sc provides the maximum 'directivity' by achieving full diversity.")

if __name__ == '__main__':
    analyze_stbc_diversity()
