import numpy as np

def solve_dft_problem():
    """
    Calculates the 4-point DFT for the specified interleaved sequence.
    """
    # Given 4-point DFTs
    X = np.array([1, 1j, -1, -1j])
    H = np.array([0, 1+1j, 1, 1-1j])

    # The resulting 4-point DFT, which we will calculate
    Z = np.zeros(4, dtype=np.complex128)

    print("Thinking Steps and Formulas:")
    print("1. The 4-point DFT Z(k) of the 8-point interleaved sequence is found by first time-aliasing the sequence.")
    print("2. The aliased 4-point sequence z(n) can be related to X(k) and H(k) in the frequency domain.")
    print("-" * 30)

    # Calculate Z(0)
    Z[0] = X[0] + H[0]
    print(f"Z[0] = X[0] + H[0] = {X[0].real} + {H[0].real} = {Z[0]}")

    # Calculate Z[1]
    Z[1] = X[2] - 1j * H[2]
    print(f"Z[1] = X[2] - j*H[2] = ({X[2]}) - 1j*({H[2]}) = {Z[1]}")

    # Calculate Z(2)
    Z[2] = X[0] - H[0]
    print(f"Z[2] = X[0] - H[0] = {X[0].real} - {H[0].real} = {Z[2]}")

    # Calculate Z(3)
    Z[3] = X[2] + 1j * H[2]
    print(f"Z[3] = X[2] + j*H[2] = ({X[2]}) + 1j*({H[2]}) = {Z[3]}")
    
    print("-" * 30)
    print("Final 4-point DFT Z(k):")
    print(Z)

solve_dft_problem()