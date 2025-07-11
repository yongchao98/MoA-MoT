import numpy as np

def format_complex(c):
    """Helper function to format complex numbers nicely for printing."""
    r, i = np.round(c.real, 2), np.round(c.imag, 2)
    
    if abs(i) < 1e-9:
        return f"{r:.0f}"
    
    if abs(r) < 1e-9:
        if abs(abs(i) - 1) < 1e-9:
            return f"j" if i > 0 else f"-j"
        else:
            return f"{i:.0f}j"

    op = "+" if i > 0 else "-"
    i_abs = abs(i)

    if abs(i_abs - 1) < 1e-9:
        return f"({r:.0f}{op}j)"
    else:
        return f"({r:.0f}{op}{i_abs:.0f}j)"


def main():
    """
    Calculates the 4-point DFT of the interleaved sequence and prints the steps.
    """
    # Given 4-point DFTs
    X = np.array([1, 1j, -1, -1j], dtype=complex)
    H = np.array([0, 1 + 1j, 1, 1 - 1j], dtype=complex)

    # Let Y_final(k) be the 4-point DFT of the interleaved sequence.
    Y_final = []

    print("To find the 4-point DFT for the 8-point sequence {x(0),h(0),...,x(3),h(3)},")
    print("we use the formula: Y(k) = X(2k mod 4) + (W_4^k) * H(2k mod 4)\n")

    # Calculate for k = 0, 1, 2, 3
    for k in range(4):
        # Calculate the twiddle factor W_4^k = exp(-j*2*pi*k/4)
        w_factor = np.exp(-1j * 2 * np.pi * k / 4)
        
        # Get the index for X and H, which is (2*k mod 4)
        idx = (2 * k) % 4
        
        # Get the corresponding values from the input DFTs
        x_val = X[idx]
        h_val = H[idx]
        
        # Calculate the final DFT value for the current k
        y_val = x_val + w_factor * h_val
        Y_final.append(y_val)
        
        print(f"For k={k}:")
        print(f"Y({k}) = X({idx}) + W_4^{k} * H({idx})")
        print(f"    = {format_complex(x_val)} + {format_complex(w_factor)} * {format_complex(h_val)}")
        print(f"    = {format_complex(y_val)}\n")

    # Format the final list for display
    final_result_list = [format_complex(c).replace('(','').replace(')','') for c in Y_final]
    print(f"The final 4-point DFT is Y(k) = [{','.join(final_result_list)}]")


if __name__ == "__main__":
    main()