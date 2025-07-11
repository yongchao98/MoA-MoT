import numpy as np

def solve_dft():
    """
    Calculates the 8-point DFT of an interleaved sequence y(n) = {x(0),h(0),...,x(3),h(3)}
    given the 4-point DFTs X(k) and H(k).
    """
    # Given 4-point DFTs
    X = np.array([1, 1j, -1, -1j])
    H = np.array([0, 1 + 1j, 1, 1 - 1j])
    
    N = 8
    Y = np.zeros(N, dtype=complex)

    print("The 8-point DFT Y(k) is calculated using the decimation-in-time FFT formulas:")
    print("For k = 0, 1, 2, 3:")
    print("Y(k)   = X(k) + W_8^k * H(k)")
    print("Y(k+4) = X(k) - W_8^k * H(k)")
    print("-" * 50)

    # Loop through the first half of the points
    for k in range(N // 2):
        # Calculate the twiddle factor W_N^k = exp(-2j * pi * k / N)
        twiddle_factor = np.exp(-2j * np.pi * k / N)
        
        # Calculate the term W_8^k * H(k)
        term = twiddle_factor * H[k]
        
        # Calculate Y(k) and Y(k + N/2)
        Y[k] = X[k] + term
        Y[k + N // 2] = X[k] - term

        # Print the detailed calculation for Y(k)
        print(f"For k={k}:")
        print(f"Y({k}) = {X[k]} + ({twiddle_factor:.4f}) * ({H[k]})")
        print(f"   = {X[k]} + ({term:.4f})")
        print(f"   = {Y[k]:.4f}\n")
        
        # Print the detailed calculation for Y(k+4)
        print(f"Y({k+4}) = {X[k]} - ({twiddle_factor:.4f}) * ({H[k]})")
        print(f"   = {X[k]} - ({term:.4f})")
        print(f"   = {Y[k+N//2]:.4f}\n")
        print("-" * 50)

    print("The final 8-point DFT vector Y(k) is:")
    # Clean up near-zero values for cleaner printing
    Y[np.isclose(Y, 0)] = 0
    print(Y)
    
    # Return the final answer in the specified format
    # The string representation is created for the final output format.
    final_answer_str = f"[{Y[0].real}, {Y[1].real:.4f}{Y[1].imag:+.4f}j, {Y[2].real:.0f}{Y[2].imag:+.0f}j, {Y[3].real:.4f}{Y[3].imag:+.4f}j, {Y[4].real}, {Y[5].real:.4f}{Y[5].imag:+.4f}j, {Y[6].real:.0f}{Y[6].imag:+.0f}j, {Y[7].real:.4f}{Y[7].imag:+.4f}j]"
    return final_answer_str

# Run the solver and print the final answer in the required format.
final_answer = solve_dft()
print(f"<<<{final_answer}>>>")
