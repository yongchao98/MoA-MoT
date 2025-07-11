import numpy as np

def solve_dft():
    """
    Calculates the 8-point DFT of an interleaved sequence from the
    4-point DFTs of the original sequences.
    """
    # Given 4-point DFTs
    X = np.array([1, 1j, -1, -1j], dtype=complex)
    H = np.array([0, 1 + 1j, 1, 1 - 1j], dtype=complex)

    # Initialize the 8-point DFT result array
    N = 8
    Y = np.zeros(N, dtype=complex)

    print("The components of the 8-point DFT Y(k) are calculated as follows:")

    # Calculate Y(k) for k = 0 to 7
    for k in range(N // 2):
        # Calculate the twiddle factor W_N^k
        twiddle_factor = np.exp(-2j * np.pi * k / N)
        
        # Apply the formula
        term = twiddle_factor * H[k]
        Y[k] = X[k] + term
        Y[k + N // 2] = X[k] - term

    # Output each number in the final resulting DFT vector Y
    # The final "equation" is Y(k) = value
    for k in range(N):
        # Format the complex number for cleaner output
        real_part = np.real(Y[k])
        imag_part = np.imag(Y[k])
        
        # Clean up floating point inaccuracies for display
        if np.isclose(real_part, 0): real_part = 0.0
        if np.isclose(imag_part, 0): imag_part = 0.0

        print(f"Y({k}) = {real_part:g} + {imag_part:g}j")

# Execute the function
solve_dft()