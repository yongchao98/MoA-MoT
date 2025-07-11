import numpy as np

def solve_dft():
    """
    Calculates the 8-point DFT of an interleaved sequence based on the
    4-point DFTs of the original sequences.
    """
    # Given 4-point DFTs
    X = np.array([1, 1j, -1, -1j], dtype=complex)
    H = np.array([0, 1 + 1j, 1, 1 - 1j], dtype=complex)

    # Initialize the 8-point DFT result vector
    Y = np.zeros(8, dtype=complex)

    print("This script calculates the 8-point DFT, Y(k), for the interleaved sequence")
    print("y(m) = {x(0), h(0), x(1), h(1), x(2), h(2), x(3), h(3)}.")
    print("\nUsing the decimation-in-time formulas:")
    print("Y(k)   = X(k) + (W_8^k) * H(k)")
    print("Y(k+4) = X(k) - (W_8^k) * H(k)")
    print("--------------------------------------------------\n")

    # Loop for k from 0 to 3 to calculate all 8 points of Y(k)
    for k in range(4):
        # Calculate the 8-point twiddle factor W_8^k
        w_8_k = np.exp(-1j * 2 * np.pi * k / 8)
        
        # Calculate the term to be added/subtracted
        term = w_8_k * H[k]
        
        # Calculate Y(k) and Y(k+4)
        Y[k] = X[k] + term
        Y[k+4] = X[k] - term

        # Print the detailed calculation for Y(k)
        print(f"For k={k}:")
        print(f"Y({k}) = X({k}) + (W_8^{k}) * H({k})")
        print(f"Y({k}) = ({X[k]:.4f}) + ({w_8_k:.4f}) * ({H[k]:.4f})")
        print(f"Y({k}) = ({X[k]:.4f}) + ({term:.4f})")
        print(f"Y({k}) = {Y[k]:.4f}\n")

        # Print the detailed calculation for Y(k+4)
        print(f"Y({k+4}) = X({k}) - (W_8^{k}) * H({k})")
        print(f"Y({k+4}) = ({X[k]:.4f}) - ({w_8_k:.4f}) * ({H[k]:.4f})")
        print(f"Y({k+4}) = ({X[k]:.4f}) - ({term:.4f})")
        print(f"Y({k+4}) = {Y[k+4]:.4f}")
        print("--------------------------------------------------")

    print("\nThe final 8-point DFT vector Y(k) is:")
    # Use np.round to clean up floating point inaccuracies for printing
    Y_rounded = np.round(Y, decimals=4)
    print(Y_rounded)

solve_dft()
<<<[1.+0.j, 1.4142+1.j, -1.-1.j, -1.4142-1.j, 1.+0.j, -1.4142+1.j, -1.+1.j, 1.4142-1.j]>>>