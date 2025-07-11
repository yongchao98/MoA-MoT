import numpy as np

def main():
    """
    Calculates the 4-point DFT of an interleaved 8-point sequence
    using the frequency-domain properties of the DFT.
    """
    # Given 4-point DFTs
    X = np.array([1, 1j, -1, -1j], dtype=complex)
    H = np.array([0, 1+1j, 1, 1-1j], dtype=complex)

    print("Given DFTs:")
    print(f"X(k) = {X}")
    print(f"H(k) = {H}\n")

    # The problem is to find the 4-point DFT, Z(k), of the 8-point sequence
    # formed by interleaving x(n) and h(n). This corresponds to taking the
    # even-indexed samples of the 8-point DFT of the interleaved sequence, Y(k).
    # Z(k) = Y(2k) for k=0,1,2,3.
    # The formulas are:
    # Z(0) = Y(0) = X(0) + W_8^0 * H(0)
    # Z(1) = Y(2) = X(2) + W_8^2 * H(2)
    # Z(2) = Y(4) = X(0) - W_8^0 * H(0)
    # Z(3) = Y(6) = X(2) - W_8^2 * H(2)
    
    # Define necessary twiddle factors
    # W_8^0 = exp(0) = 1
    # W_8^2 = exp(-j*2*pi*2/8) = exp(-j*pi/2) = -j
    W8_0 = 1.0
    W8_2 = -1.0j

    print("Calculating the 4-point DFT Z(k) step-by-step:\n")

    # Calculate Z(0)
    z0 = X[0] + W8_0 * H[0]
    print(f"Z(0) = X(0) + (W_8^0) * H(0) = {X[0].real} + ({W8_0.real}) * {H[0].real} = {z0.real}")

    # Calculate Z(1)
    z1 = X[2] + W8_2 * H[2]
    print(f"Z(1) = X(2) + (W_8^2) * H(2) = {X[2].real} + ({W8_2}) * {H[2].real} = {z1}")

    # Calculate Z(2)
    z2 = X[0] - W8_0 * H[0]
    print(f"Z(2) = X(0) - (W_8^0) * H(0) = {X[0].real} - ({W8_0.real}) * {H[0].real} = {z2.real}")

    # Calculate Z(3)
    z3 = X[2] - W8_2 * H[2]
    print(f"Z(3) = X(2) - (W_8^2) * H(2) = {X[2].real} - ({W8_2}) * {H[2].real} = {z3}")

    # Combine into the final DFT vector
    Z = np.array([z0, z1, z2, z3])

    print("\nThe final 4-point DFT for the sequence is:")
    # Print the array in a readable format
    print(f"Z(k) = [{Z[0]}, {Z[1]}, {Z[2]}, {Z[3]}]")

if __name__ == "__main__":
    main()