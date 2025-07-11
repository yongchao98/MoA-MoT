import numpy as np

def format_complex(c, precision=4):
    """A helper function to format complex numbers for printing."""
    real = c.real
    imag = c.imag
    if np.isclose(imag, 0):
        return f"{real:.{precision}f}"
    if np.isclose(real, 0):
        # Handle -0.0j case for cleaner output
        if np.isclose(imag, -0.0):
            return f"{-imag:.{precision}f}j"
        return f"{imag:.{precision}f}j"
    sign = "+" if imag > 0 else "-"
    return f"({real:.{precision}f} {sign} {abs(imag):.{precision}f}j)"

def solve_dft():
    """
    Calculates the 8-point DFT of an interleaved sequence based on the
    4-point DFTs of the original sequences.
    """
    # Given 4-point DFTs
    X = np.array([1, 1j, -1, -1j], dtype=complex)
    H = np.array([0, 1 + 1j, 1, 1 - 1j], dtype=complex)

    N = 8  # The size of the final DFT
    Y = np.zeros(N, dtype=complex)

    print("Calculating the 8-point DFT Y(k) using the formula:")
    print("Y(k)   = X(k) + W_8^k * H(k)")
    print("Y(k+4) = X(k) - W_8^k * H(k)\n")

    # Loop for k from 0 to N/2 - 1 (which is 3)
    for k in range(N // 2):
        # Calculate the twiddle factor W_N^k
        twiddle_factor = np.exp(-1j * 2 * np.pi * k / N)
        
        # The term to be added/subtracted
        term = twiddle_factor * H[k]
        
        # Calculate Y(k) and Y(k + N/2)
        Y[k] = X[k] + term
        Y[k + N//2] = X[k] - term

        print(f"--- k = {k} ---")
        print(f"W_8^{k} = {format_complex(twiddle_factor)}")
        print(f"Y({k}) = X({k}) + W_8^{k}*H({k}) = {format_complex(X[k])} + {format_complex(twiddle_factor)} * {format_complex(H[k])} = {format_complex(Y[k])}")
        print(f"Y({k+4}) = X({k}) - W_8^{k}*H({k}) = {format_complex(X[k])} - {format_complex(twiddle_factor)} * {format_complex(H[k])} = {format_complex(Y[k+N//2])}\n")

    print("--- Final Result ---")
    print("The 8-point DFT for the sequence {x(0),h(0),x(1),h(1),...} is Y(k):")
    
    # Create a nicely formatted string for the final array
    final_Y_str = "[" + ", ".join([format_complex(y) for y in Y]) + "]"
    print(final_Y_str)
    
    # This part is for the final answer block as per instructions
    # Note: The '<<<...>>>' block should be the very last thing in the response.
    # The print statement below generates the content for that block.
    # print(f"\n<<<{final_Y_str}>>>")

solve_dft()
