import numpy as np

def solve_dft():
    """
    Calculates the 8-point DFT for an interleaved sequence using the
    given 4-point DFTs of the original sequences.
    """
    # Given 4-point DFTs
    X = np.array([1, 1j, -1, -1j])
    H = np.array([0, 1+1j, 1, 1-1j])

    N = 8  # Length of the final DFT
    M = 4  # Length of the original DFTs

    # Initialize the 8-point DFT Y(k)
    Y = np.zeros(N, dtype=complex)

    print("This script calculates the 8-point DFT Y(k) for the interleaved sequence.")
    print("y(n) = {x(0), h(0), x(1), h(1), x(2), h(2), x(3), h(3)}")
    print("\nWe use the following formulas based on the decimation-in-time FFT property:")
    print("  Y(k)     = X(k) + W_8^k * H(k)")
    print("  Y(k+4)   = X(k) - W_8^k * H(k)")
    print("for k = 0, 1, 2, 3, where W_8 = exp(-2j*pi/8).\n")

    def format_complex(c):
        """Helper function to format complex numbers for clean printing."""
        real, imag = c.real, c.imag
        if np.isclose(real, 0) and np.isclose(imag, 0): return "0"
        if np.isclose(real, 0): return f"{imag:.4f}j"
        if np.isclose(imag, 0): return f"{real:.4f}"
        return f"({real:.4f} {'+' if imag >= 0 else '-'} {abs(imag):.4f}j)"

    for k in range(M):
        # Calculate the twiddle factor W_N^k
        W_N_k = np.exp(-2j * np.pi * k / N)
        
        # Calculate the product term for the equations
        term = W_N_k * H[k]
        
        # Calculate Y(k) and Y(k+4)
        Y[k] = X[k] + term
        Y[k+M] = X[k] - term

        # Print the equations with the specific numbers for this step
        print(f"--- For k = {k} ---")
        print(f"  Y({k}) = X({k}) + (W_8^{k} * H({k}))")
        print(f"       = {format_complex(X[k])} + ({format_complex(W_N_k)} * {format_complex(H[k])})")
        print(f"       = {format_complex(X[k])} + {format_complex(term)}")
        print(f"       = {format_complex(Y[k])}\n")

        print(f"  Y({k+M}) = X({k}) - (W_8^{k} * H({k}))")
        print(f"         = {format_complex(X[k])} - ({format_complex(W_N_k)} * {format_complex(H[k])})")
        print(f"         = {format_complex(X[k])} - {format_complex(term)}")
        print(f"         = {format_complex(Y[k+M])}\n")

    print("--- Final Result ---")
    print("The complete 8-point DFT vector Y(k) is:")
    result_str = f"Y = [{', '.join([format_complex(c) for c in Y])}]"
    print(result_str)

solve_dft()
<<<[1.0000, (1.4142 + 1.0000j), (-1.0000 - 1.0000j), (-1.4142 - 1.0000j), 1.0000, (-1.4142 + 1.0000j), (-1.0000 + 1.0000j), (1.4142 - 1.0000j)]>>>