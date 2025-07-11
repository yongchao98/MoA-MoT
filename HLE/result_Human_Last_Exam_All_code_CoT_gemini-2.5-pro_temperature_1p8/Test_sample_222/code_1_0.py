import numpy as np

def solve_dft():
    """
    Calculates the 8-point DFT of an interleaved sequence.
    """
    # Given 4-point DFTs
    X = np.array([1, 1j, -1, -1j], dtype=complex)
    H = np.array([0, 1 + 1j, 1, 1 - 1j], dtype=complex)
    
    # The new sequence y(n) is 8-points long, so we calculate an 8-point DFT Y(k).
    N = 8
    Y = np.zeros(N, dtype=complex)

    print("The 8-point DFT Y(k) can be found using the 4-point DFTs X(k) and H(k).")
    print("The formulas are derived from the decimation-in-time FFT algorithm:")
    print("For k = 0, 1, 2, 3: Y(k) = X(k) + W_8^k * H(k)")
    print("For k = 0, 1, 2, 3: Y(k+4) = X(k) - W_8^k * H(k)")
    print("-" * 50)
    print("Starting calculations:\n")
    
    # Loop for k from 0 to 3
    for k in range(N // 2):
        # Calculate the twiddle factor W_N^k
        W_Nk = np.exp(-2j * np.pi * k / N)
        
        # Calculate the butterfly operation
        term = W_Nk * H[k]
        Y[k] = X[k] + term
        Y[k + N // 2] = X[k] - term

        # Print the detailed calculation for Y(k)
        print(f"For k = {k}:")
        print(f"  Y({k}) = X({k}) + W_8^{k} * H({k})")
        print(f"       = {X[k]:.4f} + ({W_Nk:.4f}) * ({H[k]:.4f})")
        print(f"       = {X[k]:.4f} + ({term:.4f})")
        print(f"       = {Y[k]:.4f}\n")
        
        # Print the detailed calculation for Y(k+4)
        print(f"For k = {k+N//2}:")
        print(f"  Y({k+N//2}) = X({k}) - W_8^{k} * H({k})")
        print(f"       = {X[k]:.4f} - ({W_Nk:.4f}) * ({H[k]:.4f})")
        print(f"       = {X[k]:.4f} - ({term:.4f})")
        print(f"       = {Y[k+N//2]:.4f}\n")
        print("-" * 50)

    print("Final 8-point DFT sequence Y(k) is:")
    # Using a custom print to make it cleaner
    y_str = "[" + ", ".join([f"{val:.4f}".replace("+0.0000j", "").replace("-0.0000j", "") for val in Y]) + "]"
    print(y_str)
    
    # Return the final result in the specified format
    # np.round is used to clean up floating point inaccuracies for the final answer
    final_answer = np.round(Y, 4)
    return final_answer

# Execute the function and get the final answer string
final_Y = solve_dft()

# Format the final answer for the <<<>>> block
answer_string = "[" + ", ".join([f"{val:.4f}".replace('j', 'j') for val in final_Y]) + "]"
print(f"\n<<<{answer_string}>>>")