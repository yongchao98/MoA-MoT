import numpy as np

def solve_dft():
    """
    Calculates the 8-point DFT of an interleaved sequence given the 4-point DFTs
    of the original sequences.
    """
    # Given 4-point DFTs
    # X(k) = [1, j, -1, -j]
    # H(k) = [0, 1+j, 1, 1-j]
    X = np.array([1, 1j, -1, -1j])
    H = np.array([0, 1 + 1j, 1, 1 - 1j])

    # The sequence y(n) is {x(0), h(0), x(1), h(1), ..., x(3), h(3)}.
    # The 8-point DFT Y(k) is related to X(k) and H(k) by the formula:
    # Y(k) = X(k mod 4) + W_8^k * H(k mod 4)
    # where W_8 = exp(-j*2*pi/8)

    print("The 8-point DFT Y(k) is calculated using the formula: Y(k) = X(k mod 4) + W_8^k * H(k mod 4)\n")

    Y_result = []
    # Loop for k from 0 to 7
    for k in range(8):
        # Twiddle factor W_8^k
        W = np.exp(-1j * 2 * np.pi * k / 8)
        
        # Get the corresponding values from X and H, using modulo for periodicity
        x_val = X[k % 4]
        h_val = H[k % 4]
        
        # Calculate Y(k)
        y_val = x_val + W * h_val
        Y_result.append(y_val)
        
        # Formatting for clear output
        # Using np.round to clean up floating point inaccuracies
        x_str = str(np.round(x_val, 4)).replace('j', 'j')
        h_str = str(np.round(h_val, 4)).replace('j', 'j')
        W_str = str(np.round(W, 4)).replace('j', 'j')
        y_str = str(np.round(y_val, 4)).replace('j', 'j')

        print(f"For k = {k}:")
        print(f"Y({k}) = X({k % 4}) + W_8^{k} * H({k % 4})")
        print(f"   = {x_str} + {W_str} * {h_str}")
        print(f"   = {y_str}\n")

    # This part is for the final answer block and will not be printed by the script
    final_answer_str = ", ".join([str(np.round(y, 4)).replace('j', 'j') for y in Y_result])
    return f"[{final_answer_str}]"


# Execute the function to print the detailed calculations
if __name__ == "__main__":
    solve_dft()
