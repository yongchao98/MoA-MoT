import numpy as np

def format_complex(c):
    """Helper function to format complex numbers for clean output."""
    c = np.round(c, 10) # Clean up potential floating point inaccuracies
    if c.imag == 0:
        return f"{c.real:.4g}"
    if c.real == 0:
        if c.imag == 1:
            return "j"
        if c.imag == -1:
            return "-j"
        return f"{c.imag:.4g}j"
    sign = "+" if c.imag > 0 else "-"
    return f"({c.real:.4g} {sign} {abs(c.imag):.4g}j)"

def solve_dft():
    """
    Solves the problem of finding the 4-point DFT related to the interleaved sequence.
    """
    # Given 4-point DFTs
    X = np.array([1, 1j, -1, -1j])
    H = np.array([0, 1+1j, 1, 1-1j])

    print("Step 1: Define the input 4-point DFTs.")
    print(f"X(k) = [1, j, -1, -j]")
    print(f"H(k) = [0, 1+j, 1, 1-j]\n")
    
    print("Step 2: The new 8-point sequence y(n) is an interleaving of x(n) and h(n).")
    print("We can find its 8-point DFT, Y(k), using the Decimation-In-Time FFT formulas:")
    print("Y(k)     = X(k) + W_8^k * H(k)  (for k=0,1,2,3)")
    print("Y(k+4)   = X(k) - W_8^k * H(k)  (for k=0,1,2,3)\n")

    print("Step 3: The question 'Find 4-point DFT for...' is interpreted as finding the aliased sequence {Y(0), Y(2), Y(4), Y(6)}.\n")

    print("Step 4: Calculate the required twiddle factors W_8^k = exp(-j*2*pi*k/8).")
    # Twiddle factors
    W8_0 = np.exp(-1j * 2 * np.pi * 0 / 8)
    W8_2 = np.exp(-1j * 2 * np.pi * 2 / 8)
    print(f"W_8^0 = {format_complex(W8_0)}")
    print(f"W_8^2 = exp(-j*pi/2) = {format_complex(W8_2)}\n")

    print("Step 5: Compute each term of the resulting 4-point sequence.\n")

    # Calculate Y(0)
    Y0 = X[0] + W8_0 * H[0]
    print(f"Y(0) = X(0) + W_8^0 * H(0)")
    print(f"     = {format_complex(X[0])} + {format_complex(W8_0)} * {format_complex(H[0])}")
    print(f"     = {format_complex(X[0])} + {format_complex(W8_0 * H[0])}")
    print(f"     = {format_complex(Y0)}\n")

    # Calculate Y(2)
    Y2 = X[2] + W8_2 * H[2]
    print(f"Y(2) = X(2) + W_8^2 * H(2)")
    print(f"     = {format_complex(X[2])} + {format_complex(W8_2)} * {format_complex(H[2])}")
    print(f"     = {format_complex(X[2])} + {format_complex(W8_2 * H[2])}")
    print(f"     = {format_complex(Y2)}\n")
    
    # Calculate Y(4)
    Y4 = X[0] - W8_0 * H[0]
    print(f"Y(4) = X(0) - W_8^0 * H(0)")
    print(f"     = {format_complex(X[0])} - {format_complex(W8_0)} * {format_complex(H[0])}")
    print(f"     = {format_complex(X[0])} - {format_complex(W8_0 * H[0])}")
    print(f"     = {format_complex(Y4)}\n")

    # Calculate Y(6)
    Y6 = X[2] - W8_2 * H[2]
    print(f"Y(6) = X(2) - W_8^2 * H(2)")
    print(f"     = {format_complex(X[2])} - ({format_complex(W8_2)}) * {format_complex(H[2])}")
    print(f"     = {format_complex(X[2])} - ({format_complex(W8_2 * H[2])})")
    print(f"     = {format_complex(Y6)}\n")
    
    final_result = [Y0, Y2, Y4, Y6]
    print("The final 4-point sequence is [Y(0), Y(2), Y(4), Y(6)]")
    print("Result: [" + ", ".join([format_complex(c) for c in final_result]) + "]")

if __name__ == '__main__':
    solve_dft()
    <<<[1, (-1 - 1j), 1, (-1 + 1j)]>>>